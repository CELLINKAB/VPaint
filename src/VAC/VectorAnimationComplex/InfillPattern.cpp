#include "InfillPattern.h"

#include <optional>
#include <QDebug>
#include <QtOpenGL>

InfillPattern::InfillPattern() = default;

InfillPattern::~InfillPattern() = default;

QPointF rotateAndNormalizePoint(double angle, const QPointF& vector)
{
    constexpr auto r = 1;
    auto phi = atan2(vector.y(), vector.x());
    phi += angle;
    return QPointF{r * cos(phi), r * sin(phi)};
}

std::optional<QPointF> lineLineIntersection(QPointF A, QPointF B, QPointF C, QPointF D)
{
    // Line AB represented as a1x + b1y = c1
    const double a1 = B.y() - A.y();
    const double b1 = A.x() - B.x();
    const double c1 = a1*(A.x()) + b1*(A.y());

    // Line CD represented as a2x + b2y = c2
    const double a2 = D.y() - C.y();
    const double b2 = C.x() - D.x();
    const double c2 = a2*(C.x())+ b2*(C.y());

    const double determinant = a1*b2 - a2*b1;

    if (determinant == 0)
    {
        // The lines are parallel. This is simplified
        // by returning a pair of FLT_MAX
        return {};
    }
    else
    {
        const double x = (b2*c1 - b1*c2)/determinant;
        const double y = (a1*c2 - a2*c1)/determinant;
        return {QPointF(x, y)};
    }
}

bool ccw(QPointF A, QPointF B, QPointF C)
{
    return (C.y() - A.y()) * (B.x() - A.x()) > (B.y() - A.y()) * (C.x() - A.x());
}

bool intersect(QPointF A, QPointF B, QPointF C, QPointF D)
{
    return ccw(A, C, D) != ccw(B, C, D) && ccw(A, B, C) != ccw(A, B, D);
}

bool isClockwise(const QPolygonF& polygon)
{
    auto sum = 0;
    for (int i = 0; i < polygon.size(); i++) {
        const auto firstPoint = polygon[i];
        const auto secondPoint = polygon[(i + 1) % polygon.size()];
        sum += (secondPoint.x() - firstPoint.x()) * (firstPoint.y() + secondPoint.y());
    }
    return sum < 0;
}

std::optional<QPointF> segmentSegmentIntersection(QPointF A, QPointF B, QPointF C, QPointF D)
{
    if (intersect(A, B, C, D)
            || pointToLineSegmentDistanceSquared(A, C, D) == 0
            || pointToLineSegmentDistanceSquared(B, C, D) == 0
            || pointToLineSegmentDistanceSquared(C, A, B) == 0
            || pointToLineSegmentDistanceSquared(D, A, B) == 0) {
        const auto intersection = lineLineIntersection(A, B, C, D);
        return intersection;
    }
    return {};
}

QVector<QPair<QPointF, int>> linePolygonIntersections(QPointF A, QPointF B, QPolygonF inset)
{
    QVector<QPair<QPointF, int>> intersectedLines{};
    for (int i = 0; i < inset.size(); i++) {
        if (const auto point = segmentSegmentIntersection(A, B, inset[i], inset[(i + 1) % inset.size()])) {
            intersectedLines << QPair{point.value(), i};
        }
    }
    return intersectedLines;
}

bool QPointFLessThan(const QPair<QPointF, int> &p1, const QPair<QPointF, int> &p2)
{
    return std::pow(p1.first.x(), 2) + std::pow(p1.first.y(), 2) < std::pow(p2.first.x(), 2) + std::pow(p2.first.y(), 2);
}

QVector<QPair<QPointF, int>> sortPointsByDistance(QPointF startPosition, const QVector<QPair<QPointF, int>>& foundIntersectionPoints)
{
    QVector<QPair<QPointF, int>> subtractedIntersectionPoints;
    for (const auto foundIntersectionPoint: foundIntersectionPoints) {
        subtractedIntersectionPoints << QPair{foundIntersectionPoint.first - startPosition, foundIntersectionPoint.second};
    }
    std::sort(subtractedIntersectionPoints.begin(), subtractedIntersectionPoints.end(), QPointFLessThan);
    QVector<QPair<QPointF, int>> addedBackIntersectionPoints;
    for (const auto& intersection: subtractedIntersectionPoints) {
        addedBackIntersectionPoints << QPair{intersection.first + startPosition, intersection.second};
    }
    return addedBackIntersectionPoints;
}

double pointToLineSegmentDistanceSquared(QPointF p, QPointF startPoint, QPointF endPoint)
{
    const auto x = p.x();
    const auto y = p.y();
    const auto x1 = startPoint.x();
    const auto y1 = startPoint.y();
    const auto x2 = endPoint.x();
    const auto y2 = endPoint.y();

    const auto A = x - x1;
    const auto B = y - y1;
    const auto C = x2 - x1;
    const auto D = y2 - y1;

    const auto dot = A * C + B * D;
    const auto lengthSquared = C * C + D * D;
    auto param = -1.0;
    if (lengthSquared != 0) {
        param = dot / lengthSquared;
    }

    auto xx = x1 + param * C;
    auto yy = y1 + param * D;
    if (param < 0.0) {
        xx = x1;
        yy = y1;
    } else if (param > 1.0) {
        xx = x2;
        yy = y2;
    }
    const auto dx = x - xx;
    const auto dy = y - yy;
    return std::pow(dx, 2) + std::pow(dy, 2);
}

double polygonCircumferenceDistance(const QPolygonF& polygon)
{
    auto distance = 0.0;
    for (int i = 0; i < polygon.size() - 1; i++) {
        const auto firstPoint = polygon[i];
        const auto secondPoint = polygon[i + 1];
        distance += std::sqrt(std::pow(secondPoint.x() - firstPoint.x(), 2) + std::pow(secondPoint.y() - firstPoint.y(), 2));
    }
    return distance;
}

double pointToPointDistanceSquared(QPointF a, QPointF b) {
    return std::pow(a.x() - b.x(), 2) + std::pow(a.y() - b.y(), 2);
}

bool compare(double d1, double d2, quint8 precision)
{
    return std::abs(d1 - d2) < std::pow(10, -precision);
}

int roundDouble(double value)
{
    return (int) (value * 1000);
}

QPolygonF traverseFromStartToEnd(QPointF startPoint, QPointF endPoint, int insetStartIndex, int insetEndIndex, const QPolygonF& inset)
{
    QPolygonF path{};
    path << startPoint;
    // If we are moving in direction from start to end add endPoint
    if (insetStartIndex == insetEndIndex &&
        pointToPointDistanceSquared(startPoint, inset[insetStartIndex]) < pointToPointDistanceSquared(endPoint, inset[insetStartIndex])) {
        path<< endPoint;
    } else {
        auto insetIndex = (insetStartIndex + 1) % inset.size();
        while (insetIndex != insetStartIndex) {
            const auto startInsetPoint = inset[insetIndex];
            const auto endInsetPoint = inset[(insetIndex + 1) % inset.size()];
            if (compare(pointToLineSegmentDistanceSquared(endPoint, startInsetPoint, endInsetPoint), 0.0, 5)) {
                if (path.last() != startInsetPoint) {
                    path << startInsetPoint;
                }
                if (path.last() != endPoint) {
                    path << endPoint;
                }
                break;
            } else {
                if (path.last() != startInsetPoint) {
                    path << startInsetPoint;
                }
                if (path.last() != endInsetPoint) {
                    path << endInsetPoint;
                }
            }
            insetIndex = (insetIndex + 1) % inset.size();
        }
    }
    return path;
}

QPolygonF reversePolygonFOrientation(const QPolygonF& polygon)
{
    QPolygonF reversedPolygon{};
    for (int i = polygon.size() - 1; i >= 0; i--) {
        reversedPolygon << polygon[i];
    }

    return reversedPolygon;
}

std::tuple<int, int> startAndEndIndex(const QPolygonF& polygon, QPointF startPoint, QPointF endPoint)
{
    auto startIndex = -1;
    auto endIndex = -1;

    for (int i = 0; i < polygon.size(); i++) {
        const auto lineStartPoint = polygon[i];
        const auto lineEndPoint = polygon[(i + 1) % polygon.size()];
        if (startPoint != lineEndPoint && compare(pointToLineSegmentDistanceSquared(startPoint, lineStartPoint, lineEndPoint), 0.0, 5)) {
            startIndex = i;
        }
        if (endPoint != lineEndPoint && compare(pointToLineSegmentDistanceSquared(endPoint, lineStartPoint, lineEndPoint), 0.0, 5)) {
            endIndex = i;
        }
    }
    return std::tuple{startIndex, endIndex};
}

QVector<QVector<QPair<QPointF, int>>> pruneInfill(const QPolygonF& naiveInfill, const QPolygonF& inset)
{
    QVector<QVector<QPair<QPointF, int>>> prunedInfill{};
    QVector<QPair<QPointF, int>> connectedPrunedInfill{};
    std::optional<QPointF> lastIntersectionPoint{};
    for (int i = 0; i < naiveInfill.size() - 1; i++) {
        const auto firstPoint = naiveInfill[i];
        const auto secondPoint = naiveInfill[i + 1];

        auto intersections = linePolygonIntersections(firstPoint, secondPoint, inset);

        // Make sure to not step out and in on same interection point
        if (!intersections.isEmpty() && lastIntersectionPoint && intersections.first().first == lastIntersectionPoint.value()) {
            intersections.pop_front();
        }

        // We first have to step inside
        if (connectedPrunedInfill.isEmpty() && intersections.isEmpty()) {
            continue;
        }
        // Since we are inside and have no intersections, both points are inside
        if (intersections.isEmpty()) {
            if (firstPoint != connectedPrunedInfill.last().first) {
                connectedPrunedInfill << QPair(firstPoint, -1);
            }
            if (secondPoint != connectedPrunedInfill.last().first) {
                connectedPrunedInfill << QPair(secondPoint, -1);
            }
        } else {
            const auto sortedIntersections = sortPointsByDistance(firstPoint, intersections);
            // Stepping inside
            if (connectedPrunedInfill.isEmpty()) {
                for (int j = 0; j < sortedIntersections.size(); j++) {
                    connectedPrunedInfill << sortedIntersections[j];
                    if (j % 2 == 1) {
                        prunedInfill << connectedPrunedInfill;
                        connectedPrunedInfill.clear();
                    }
                }

            }
            // Stepping outside
            else {
                if (firstPoint != connectedPrunedInfill.last().first) {
                    connectedPrunedInfill << QPair(firstPoint, -1);
                }
                for (int j = 0; j < sortedIntersections.size(); j++) {
                    connectedPrunedInfill << sortedIntersections[j];
                    if (j % 2 == 0) {
                        prunedInfill << connectedPrunedInfill;
                        connectedPrunedInfill.clear();
                    }
                }
            }
            lastIntersectionPoint = {sortedIntersections.last().first};
        }
    }
    return prunedInfill;
}

QVector<QPointF> connectInfillAlongInset(const QVector<QVector<QPair<QPointF, int>>>& prunedInfill, const QPolygonF& inset, bool connectBackToBeginning)
{
    QVector<QPointF> connectedInfill{};
    if (!prunedInfill.empty()) {
        for (const auto pI: prunedInfill[0]) {
            connectedInfill << pI.first;
        }
    }
    auto subtractionStopFromEnd = 1;
    if (connectBackToBeginning) {
        subtractionStopFromEnd = 0;
    }
    const auto reversedPolygon = reversePolygonFOrientation(inset);
    for (int i = 0; i < prunedInfill.size() - subtractionStopFromEnd; i++) {
        const auto startConnectedInfill = prunedInfill[i];
        const auto endConnectedInfill = prunedInfill[(i + 1) % prunedInfill.size()];
        const auto startPointPair = startConnectedInfill.last();
        const auto endPointPair = endConnectedInfill.first();
        const auto startPoint = startPointPair.first;
        const auto endPoint = endPointPair.first;
        auto insetStartIndex = startPointPair.second;
        auto insetEndIndex = endPointPair.second;

        auto possiblePath = traverseFromStartToEnd(startPoint, endPoint, insetStartIndex, insetEndIndex, inset);
        const auto startAndEndIndexTuple = startAndEndIndex(reversedPolygon, startPoint, endPoint);
        insetStartIndex = std::get<0>(startAndEndIndexTuple);
        insetEndIndex = std::get<1>(startAndEndIndexTuple);
        auto possiblePathReverse = traverseFromStartToEnd(startPoint, endPoint, insetStartIndex, insetEndIndex, reversedPolygon);

        if (polygonCircumferenceDistance(possiblePath) <= polygonCircumferenceDistance(possiblePathReverse)) {
            if (connectedInfill.last() == possiblePath.first()) {
                possiblePath.pop_front();
            }
            connectedInfill << possiblePath;
        } else {
            if (connectedInfill.last() == possiblePathReverse.first()) {
                possiblePathReverse.pop_front();
            }
            connectedInfill << possiblePathReverse;
        }
        for (const auto pI: endConnectedInfill) {
            if (pI.first != connectedInfill.last()) {
                connectedInfill << pI.first;
            }
        }
    }
    return connectedInfill;
}

bool previousBorderContainsNewBorder(const QPolygonF& newBorder, const QPolygonF& previousBorder, int smallestDistanceAllowedSquared)
{
    // Make distance slightly smaller to allow points on edge
    smallestDistanceAllowedSquared *= 0.99999;
    for (const auto point : newBorder) {
        for (int i = 0; i < previousBorder.size(); i++) {
            const auto startPoint = previousBorder[i];
            const auto endPoint = previousBorder[(i + 1) % previousBorder.size()];
            if (pointToLineSegmentDistanceSquared(point, startPoint, endPoint) < smallestDistanceAllowedSquared) {
                return false;
            }
        }
    }
    return true;
}

QPointF getConvexPolygonCenterOfMass(const QPolygonF &polygon){
    QPointF centerPoint;
    for (const auto &point : polygon) {
        centerPoint += point;
    }
    centerPoint /= polygon.size();
    return centerPoint;
}

QPolygonF removeSelfIntersections(const QPolygonF& polygon)
{
    QPolygonF polygonWithoutIntersections;

    // Take a point inside the polygon
    const auto interiorPoint = getConvexPolygonCenterOfMass(polygon);
    constexpr auto maxNumber = std::numeric_limits<double>::max();

    // Intersect interior point to canvas edges with all segments and store intersections
    QVector<QPointF> farAwayPoints{{0,0}, {0, maxNumber}, {maxNumber, 0}, {maxNumber, maxNumber}};
    QVector<std::tuple<QPointF, int>> intersections;
    for (const auto farAwayPoint : farAwayPoints) {
        for (int i = 0; i < polygon.size(); i++) {
            const auto startPoint = polygon[i];
            const auto endPoint = polygon[(i + 1) % polygon.size()];
            if (const auto intersectionPoint = segmentSegmentIntersection(interiorPoint, farAwayPoint, startPoint, endPoint)){
                intersections.append({intersectionPoint.value(), i});
            }
        }
    }

    // Find the closest intersection segment
    // This is to make sure we don't end up on a segment that should be removed
    std::optional<std::tuple<QPointF, int>> closestSegmentToInteriorPoint{};
    auto smallestDistance = maxNumber;
    for (const auto &point : intersections){
        if (const auto distanceSquared = pointToPointDistanceSquared(interiorPoint, std::get<0>(point)) < smallestDistance){
            smallestDistance = distanceSquared;
            closestSegmentToInteriorPoint = point;
        }
    }
    if (!closestSegmentToInteriorPoint) {
        qWarning() << "Failed to find closestSegmentToInteriorPoint";
        return polygon;
    }

    // Start from closest segment and travel across polygon, piecing together new polygon without intersections
    auto segmentIndex = std::get<1>(closestSegmentToInteriorPoint.value());
    const auto endSegmentIndex = segmentIndex;
    const auto polygonSize = polygon.size();

    auto startPoint = std::get<0>(closestSegmentToInteriorPoint.value());// polygon[segmentIndex % polygonSize];
    do {
        // Check all lines against this line to find intersections
        const auto endPoint = polygon[(segmentIndex + 1) % polygonSize];
        const auto compareDirectionX = roundDouble(endPoint.x()) - roundDouble(startPoint.x()) > 0;
        const auto compareDirectionY = roundDouble(endPoint.y()) - roundDouble(startPoint.y()) > 0;
        auto innerSegmentIndex = (segmentIndex + 1) % polygonSize;
        // Find self intersections
        QVector<std::tuple<QPointF, int>> selfIntersections;
        do {
            if (const auto intersection = segmentSegmentIntersection(startPoint,
                                                                     endPoint,
                                                                     polygon[innerSegmentIndex % polygonSize],
                                                                     polygon[(innerSegmentIndex + 1) % polygonSize])) {
                // Make sure the point is in right direction
                const auto directionX = roundDouble(intersection->x()) - roundDouble(startPoint.x()) > 0;
                const auto directionY = roundDouble(intersection->y()) - roundDouble(startPoint.y()) > 0;
                if (intersection.value() != startPoint && directionX == compareDirectionX && directionY == compareDirectionY) {
                    selfIntersections.append({intersection.value(), innerSegmentIndex});
                }
            }
            innerSegmentIndex = (innerSegmentIndex + 1) % polygonSize;
        } while(innerSegmentIndex != segmentIndex);

        // Find closest intersection segment
        std::optional<std::tuple<QPointF, int>> closestSegmentToSegment{};
        auto smallestDistance = maxNumber;
        for (const auto &point : selfIntersections){
            if (const auto distanceSquared = pointToPointDistanceSquared(startPoint, std::get<0>(point)); distanceSquared < smallestDistance){
                smallestDistance = distanceSquared;
                closestSegmentToSegment = point;
            }
        }

        // Append correct point and increase segmentIndex
        if (!closestSegmentToSegment) {
            break;
        }
        const auto intersectionPoint = std::get<0>(closestSegmentToSegment.value());
        polygonWithoutIntersections << intersectionPoint;
        segmentIndex = std::get<1>(closestSegmentToSegment.value());
        startPoint = intersectionPoint;
    } while (segmentIndex != endSegmentIndex);

    return polygonWithoutIntersections;
}

void InfillPattern::rectiLinearVerticalInfill()
{
    QVector<QPointF> infillPoints{};

    const auto outerBoundingBox = inputPolygon_.boundingRect();
    const auto startBoundary = outerBoundingBox.left();
    const auto endBoundary = outerBoundingBox.right();
    const auto top = outerBoundingBox.top();
    const auto bottom = outerBoundingBox.bottom();
    const auto boundingPolygon = QPolygonF{{outerBoundingBox.topLeft(), outerBoundingBox.topRight(), outerBoundingBox.bottomRight(), outerBoundingBox.bottomLeft()}};

    auto directionUp = true;

    for (auto x = startBoundary; x <= endBoundary; x += spacing_) {
        auto potentialTopPoint = QPointF{x, top};
        while((!boundingPolygon.containsPoint(potentialTopPoint, Qt::OddEvenFill)) && potentialTopPoint.y() <= bottom) {
            potentialTopPoint += QPointF{0, 1};
        }
        auto potentialBottomPoint = QPointF{x, bottom};
        while((!boundingPolygon.containsPoint(potentialBottomPoint, Qt::OddEvenFill)) && potentialBottomPoint.y() >= top) {
            potentialBottomPoint -= QPointF{0, 1};
        }
        if (directionUp) {
            infillPoints << potentialBottomPoint;
            infillPoints << potentialTopPoint;
        } else {
            infillPoints << potentialTopPoint;
            infillPoints << potentialBottomPoint;
        }
        directionUp = !directionUp;
    }

    const auto prunedInfill = pruneInfill(infillPoints, inset_);

    data_ << connectInfillAlongInset(prunedInfill, inset_, false);
}

void InfillPattern::rectiLinearHorizontalInfill()
{
    QVector<QPointF> infillPoints{};

    const auto outerBoundingBox = inputPolygon_.boundingRect();
    const auto startBoundary = outerBoundingBox.top();
    const auto endBoundary = outerBoundingBox.bottom();
    const auto left = outerBoundingBox.left();
    const auto right = outerBoundingBox.right();
    const auto boundingPolygon = QPolygonF{{outerBoundingBox.topLeft(), outerBoundingBox.topRight(), outerBoundingBox.bottomRight(), outerBoundingBox.bottomLeft()}};

    auto directionRight = true;

    for (auto y = startBoundary; y <= endBoundary; y += spacing_) {
        auto potentialLeftPoint = QPointF{left, y};
        while((!boundingPolygon.containsPoint(potentialLeftPoint, Qt::OddEvenFill)) && potentialLeftPoint.x() <= right) {
            potentialLeftPoint += QPointF{1, 0};
        }
        auto potentialRightPoint = QPointF{right, y};
        while(!boundingPolygon.containsPoint(potentialRightPoint, Qt::OddEvenFill) && potentialRightPoint.x() >= left) {
            potentialRightPoint -= QPointF{1, 0};
        }
        if (directionRight) {
            infillPoints << potentialRightPoint;
            infillPoints << potentialLeftPoint;
        } else {
            infillPoints << potentialLeftPoint;
            infillPoints << potentialRightPoint;
        }
        directionRight = !directionRight;
    }

    const auto prunedInfill = pruneInfill(infillPoints, inset_);

    data_ << connectInfillAlongInset(prunedInfill, inset_, false);
}

void InfillPattern::gridInfill()
{
    rectiLinearHorizontalInfill();
    const auto outerBoundingBox = inputPolygon_.boundingRect();
    const auto startBoundary = outerBoundingBox.left();
    const auto endBoundary = outerBoundingBox.right();
    const auto top = outerBoundingBox.top();
    const auto bottom = outerBoundingBox.bottom();
    const auto boundingPolygon = QPolygonF{{outerBoundingBox.topLeft(), outerBoundingBox.topRight(), outerBoundingBox.bottomRight(), outerBoundingBox.bottomLeft()}};

    for (auto x = startBoundary; x <= endBoundary; x += spacing_) {
        auto potentialTopPoint = QPointF{x, top};
        while((!boundingPolygon.containsPoint(potentialTopPoint, Qt::OddEvenFill)) && potentialTopPoint.y() <= bottom) {
            potentialTopPoint += QPointF{0, 1};
        }
        auto potentialBottomPoint = QPointF{x, bottom};
        while((!boundingPolygon.containsPoint(potentialBottomPoint, Qt::OddEvenFill)) && potentialBottomPoint.y() >= top) {
            potentialBottomPoint -= QPointF{0, 1};
        }

        pruneAndAddLine(potentialTopPoint, potentialBottomPoint);
    }
}

void InfillPattern::concentricInfill()
{
    const auto insetSpacing = spacing_ / 2;
    const auto insetSpacingSquared = std::pow(insetSpacing, 2);
    auto currentBorder = inset_;

    // Temporary fix to prevent the endless loop which has appeared sometimes and caused the freezing and crash
    // when the concentric infill has applied and almost always on resize the shape with infill
    // TODO: Need to improve this logics because the concentric infil can be incorrect behaviour for not symmetric shapes
    // It is most noticeable on ellipses, for the right circle it's working fine
    const int max_data_count = currentBorder.boundingRect().width() < currentBorder.boundingRect().height() ?
                               currentBorder.boundingRect().width() / spacing_:
                               currentBorder.boundingRect().height() / spacing_;

    while((currentBorder.boundingRect().right() - currentBorder.boundingRect().left()) > insetSpacing &&
          (currentBorder.boundingRect().bottom() - currentBorder.boundingRect().top()) > insetSpacing &&
           data_.count() < max_data_count) {
        currentBorder << currentBorder.first();
        data_ << currentBorder;
        currentBorder.pop_back();
        const auto previousBorder = currentBorder;
        currentBorder = insetPolygon(currentBorder, insetSpacing);
        if (!previousBorderContainsNewBorder(currentBorder, previousBorder, insetSpacingSquared)) {
            break;
        }
    }
}

void InfillPattern::linearInfill()
{
    rectiLinearHorizontalInfill();
    const auto outerBoundingBox = inputPolygon_.boundingRect();
    const auto top = outerBoundingBox.top();
    const auto bottom = outerBoundingBox.bottom();
    const auto boundingPolygon = QPolygonF{{outerBoundingBox.topLeft(), outerBoundingBox.topRight(), outerBoundingBox.bottomRight(), outerBoundingBox.bottomLeft()}};

    const auto startPositionX = outerBoundingBox.left() + outerBoundingBox.width() / 2;
    const auto startPositionY = outerBoundingBox.top() + outerBoundingBox.height() / 2;

    auto potentialTopPoint = QPointF{startPositionX, startPositionY};
    while(boundingPolygon.containsPoint(potentialTopPoint, Qt::OddEvenFill) && potentialTopPoint.y() >= top) {
        potentialTopPoint += QPointF{0.7, -0.7};
    }
    potentialTopPoint -= QPointF{0.7, -0.7};
    auto potentialBottomPoint = QPointF{startPositionX, startPositionY};
    while(boundingPolygon.containsPoint(potentialBottomPoint, Qt::OddEvenFill) && potentialBottomPoint.y() <= bottom) {
        potentialBottomPoint += QPointF{-0.7, 0.7};
    }
    potentialBottomPoint -= QPointF{-0.7, 0.7};

    pruneAndAddLine(potentialTopPoint, potentialBottomPoint);

    for (auto distance = spacing_; distance <= outerBoundingBox.width() / 2; distance += spacing_) {
        auto potentialRightTopPoint = QPointF{startPositionX + distance, startPositionY + distance * 0.7};
        while(boundingPolygon.containsPoint(potentialRightTopPoint, Qt::OddEvenFill) && potentialRightTopPoint.y() >= top) {
            potentialRightTopPoint += QPointF{0.7, -0.7};
        }
        potentialRightTopPoint -= QPointF{0.7, -0.7};
        auto potentialRightBottomPoint = QPointF{startPositionX + distance, startPositionY + distance * 0.7};
        while(boundingPolygon.containsPoint(potentialRightBottomPoint, Qt::OddEvenFill) && potentialRightBottomPoint.y() <= bottom) {
            potentialRightBottomPoint += QPointF{-0.7, 0.7};
        }
        potentialRightBottomPoint -= QPointF{-0.7, 0.7};

        pruneAndAddLine(potentialRightTopPoint, potentialRightBottomPoint);

        auto potentialLeftTopPoint = QPointF{startPositionX - distance, startPositionY - distance * 0.7};
        while(boundingPolygon.containsPoint(potentialLeftTopPoint, Qt::OddEvenFill) && potentialLeftTopPoint.y() >= top) {
            potentialLeftTopPoint += QPointF{0.7, -0.7};
        }
        potentialLeftTopPoint -= QPointF{0.7, -0.7};
        auto potentialLeftBottomPoint = QPointF{startPositionX - distance, startPositionY - distance * 0.7};
        while(boundingPolygon.containsPoint(potentialLeftBottomPoint, Qt::OddEvenFill) && potentialLeftBottomPoint.y() <= bottom) {
            potentialLeftBottomPoint += QPointF{-0.7, 0.7};
        }
        potentialLeftBottomPoint -= QPointF{-0.7, 0.7};

        pruneAndAddLine(potentialLeftTopPoint, potentialLeftBottomPoint);
    }

    for (auto distance = spacing_ / 2; distance <= outerBoundingBox.width() / 2; distance += spacing_) {
        auto potentialRightTopPoint = QPointF{startPositionX + distance, startPositionY - distance * 0.7};
        while(boundingPolygon.containsPoint(potentialRightTopPoint, Qt::OddEvenFill) && potentialRightTopPoint.y() >= top) {
            potentialRightTopPoint += QPointF{-0.7, -0.7};
        }
        potentialRightTopPoint -= QPointF{-0.7, -0.7};
        auto potentialRightBottomPoint = QPointF{startPositionX + distance, startPositionY - distance * 0.7};
        while(boundingPolygon.containsPoint(potentialRightBottomPoint, Qt::OddEvenFill) && potentialRightBottomPoint.y() <= bottom) {
            potentialRightBottomPoint += QPointF{0.7, 0.7};
        }
        potentialRightBottomPoint -= QPointF{0.7, 0.7};

        pruneAndAddLine(potentialRightTopPoint, potentialRightBottomPoint);

        auto potentialLeftTopPoint = QPointF{startPositionX - distance, startPositionY + distance * 0.7};
        while(boundingPolygon.containsPoint(potentialLeftTopPoint, Qt::OddEvenFill) && potentialLeftTopPoint.y() >= top) {
            potentialLeftTopPoint += QPointF{-0.7, -0.7};
        }
        potentialLeftTopPoint -= QPointF{-0.7, -0.7};
        auto potentialLeftBottomPoint = QPointF{startPositionX - distance, startPositionY + distance * 0.7};
        while(boundingPolygon.containsPoint(potentialLeftBottomPoint, Qt::OddEvenFill) && potentialLeftBottomPoint.y() <= bottom) {
            potentialLeftBottomPoint += QPointF{0.7, 0.7};
        }
        potentialLeftBottomPoint -= QPointF{0.7, 0.7};

        pruneAndAddLine(potentialLeftTopPoint, potentialLeftBottomPoint);
    }
}

void InfillPattern::gyroidInfill()
{
    QVector<QPointF> infillPoints{};

    const auto outerBoundingBox = inputPolygon_.boundingRect();
    const auto startBoundary = outerBoundingBox.left();
    const auto endBoundary = outerBoundingBox.right();
    const auto top = outerBoundingBox.top();
    const auto bottom = outerBoundingBox.bottom();
    auto up = true;

    auto phaseIncrement = 1;
    for (auto x = startBoundary; x <= endBoundary; x += spacing_) {
        const auto decrementWith = std::min(0.1 * spacing_, 1.0);
        if (up) {
            for (auto y = bottom; y >= top; y -= decrementWith) {
                const auto potentialX = x + (spacing_ / 4) * sin((y - phaseIncrement * spacing_ / 4) * 4 / spacing_);
                auto potentialPoint = QPointF{potentialX, y};
                infillPoints << potentialPoint;
            }
        } else {
            for (auto y = top; y <= bottom; y += decrementWith) {
                const auto potentialX = x + (spacing_ / 4) * sin((y - phaseIncrement * spacing_ / 4) * 4 / spacing_);
                auto potentialPoint = QPointF{potentialX, y};
                infillPoints << potentialPoint;
            }
        }

        up = !up;
        phaseIncrement = phaseIncrement + 1 % 2;
    }

    const auto prunedInfill = pruneInfill(infillPoints, inset_);

    data_ << connectInfillAlongInset(prunedInfill, inset_, false);
}

void InfillPattern::honeycombInfill()
{
    QVector<QVector<QPointF>> infillPoints{};

    const auto outerBoundingBox = inputPolygon_.boundingRect();
    const auto startBoundary = outerBoundingBox.left();
    const auto endBoundary = outerBoundingBox.right();
    const auto top = outerBoundingBox.top();
    const auto bottom = outerBoundingBox.bottom();

    const auto L = spacing_;
    const auto d = 4;
    const auto rightDown = QPointF{L * std::sin(M_PI/6), L * std::cos(M_PI/6) - d / 2};
    const auto leftDown = QPointF{-L * std::sin(M_PI/6), L * std::cos(M_PI/6) - d / 2};
    const auto leftUp = QPointF{-L * std::sin(M_PI/6), -L * std::cos(M_PI/6) + d / 2};
    const auto rightUp = QPointF{L * std::sin(M_PI/6), -L * std::cos(M_PI/6) + d / 2};
    const auto boundingPolygon = QPolygonF{{outerBoundingBox.topLeft(), outerBoundingBox.topRight(), outerBoundingBox.bottomRight(), outerBoundingBox.bottomLeft()}};

    for (auto y = top; y <= bottom; y+= (rightDown.y() + d) * 2) {
        QPolygonF contiguousData;
        QPointF currentPosition{startBoundary, y};
        while (currentPosition.x() <= endBoundary) {
            contiguousData << currentPosition;
            currentPosition.setX(currentPosition.x() + L);
            contiguousData << currentPosition;
            currentPosition += rightDown;
            contiguousData << currentPosition;
            currentPosition.setX(currentPosition.x() + L);
            contiguousData << currentPosition;
            currentPosition += rightUp;
            contiguousData << currentPosition;
        }
        currentPosition = QPointF{currentPosition.x() + L, y + rightDown.y() * 2 + d};
        while (currentPosition.x() >= startBoundary) {
            contiguousData << currentPosition;
            currentPosition.setX(currentPosition.x() - L);
            contiguousData << currentPosition;
            currentPosition += leftUp;
            contiguousData << currentPosition;
            currentPosition.setX(currentPosition.x() - L);
            contiguousData << currentPosition;
            currentPosition += leftDown;
            contiguousData << currentPosition;
        }
        if (contiguousData.size() > 2) {
            contiguousData << contiguousData.first();
            infillPoints << contiguousData;
        }
    }

    for (const auto& infillPoint : infillPoints) {
        const auto prunedInfill = pruneInfill(infillPoint, inset_);
        data_ << connectInfillAlongInset(prunedInfill, inset_, true);
    }
}

void InfillPattern::pruneAndAddLine(QPointF firstPoint, QPointF secondPoint)
{
    QPolygonF line;
    line << firstPoint << secondPoint;
    const auto prunedInfill = pruneInfill(line, inset_);
    for (const auto& infill: prunedInfill) {
        QVector<QPointF> infillList{};
        for (const auto& infillPair: infill) {
            infillList << infillPair.first;
        }
        if (infillList.size() > 1) {
            data_ << infillList;
        }
    }
}

QPolygonF InfillPattern::insetPolygon(QPolygonF polygon, double distance)
{
    QPolygonF unprunedInset{};
    for (int i = 0; i < polygon.size(); i++) {
        const auto startPoint = polygon.at(i);
        const auto endPoint = polygon.at((i + 1) % (polygon.size()));
        QPointF slope{endPoint.x() - startPoint.x(), endPoint.y() - startPoint.y()};
        const auto normal = rotateAndNormalizePoint(M_PI_2, slope);
        unprunedInset << startPoint + QPointF{normal.x(), normal.y()} * distance <<  endPoint + QPointF{normal.x(), normal.y()} * distance;
    }
    for (int i = 0; i < unprunedInset.size(); i +=2) {
        // end and start
        const auto A = unprunedInset[i];
        const auto B = unprunedInset[i + 1];
        const auto C = unprunedInset[(i + 2) % unprunedInset.size()];
        const auto D = unprunedInset[(i + 3) % unprunedInset.size()];
        // If this is false we have a concave shape and didn't find the intersection
        if (const auto intersectionPoint = lineLineIntersection(A, B, C, D)) {
            unprunedInset[i + 1] = intersectionPoint.value();
            unprunedInset[(i + 2) % unprunedInset.size()] = intersectionPoint.value();
        }
    }
    QPointF currentPos{};
    QPolygonF prunedInset{};
    for (int i = 0; i < unprunedInset.size() - 1; i++) {
        const auto point = unprunedInset[i % unprunedInset.size()];
        if (point != currentPos) {
            prunedInset << point;
            currentPos = point;
        }
    }
    return removeSelfIntersections(prunedInset);
}

void InfillPattern::update(QPolygonF &polygon)
{
    spacing_ = (100 / density_) * 12;
    data_.clear();
    inset_.clear();
    QPolygonF clockwisePolygon;
    if (isClockwise(polygon)) {
        clockwisePolygon = polygon;
    } else {
        for (int i = polygon.size() - 1; i >= 0; i--) {
            clockwisePolygon << polygon[i];
        }
    }
    inset_ = insetPolygon(clockwisePolygon, insetDistance_);
    inputPolygon_ = clockwisePolygon;
    switch (pattern_) {
    case Pattern::None:
        break;
    case Pattern::Grid:
        gridInfill();
        break;
    case Pattern::RectiLinear:
        rectiLinearVerticalInfill();
        break;
    case Pattern::Concentric:
        concentricInfill();
        break;
    case Pattern::Linear:
        linearInfill();
        break;
    case Pattern::Gyroid:
        gyroidInfill();
        break;
    case Pattern::HoneyComb: {
        honeycombInfill();
        break;
    }
    default:
        qWarning() << "Unknown pattern_" << static_cast<uint8_t>(pattern_);
        break;
    }
}

const QVector<QVector<QPointF>>& InfillPattern::getPoints() const
{
    return data_;
}

void InfillPattern::draw()
{
    for (const auto& contiguousData : data_) {
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i < contiguousData.size(); i++) {
            glVertex2d(contiguousData[i].x(), contiguousData[i].y());
        }

        glEnd();
    }

    // draw inset
    #if 0
    glBegin(GL_LINE_STRIP);
    glColor3b(120, 100, 30);
    for (int i = 0; i < inset_.size(); i++) {
        glVertex2d(inset_[i].x(), inset_[i].y());
    }
    glEnd();
    #endif
}

InfillPattern::Pattern InfillPattern::pattern() const
{
    return pattern_;
}

void InfillPattern::setPattern(InfillPattern::Pattern pattern)
{
    pattern_ = pattern;
}

bool InfillPattern::isApplyInfill() const
{
    return applyInfill_;
}

void InfillPattern::setApplyInfill(bool apply)
{
    applyInfill_ = apply;
}

int InfillPattern::density() const
{
    return density_;
}

void InfillPattern::setDensity(double density)
{
    density_ = density;
}
