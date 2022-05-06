#include "InfillPattern.h"

#include <optional>

#include <QtAlgorithms>

InfillPattern::InfillPattern() = default;

InfillPattern::~InfillPattern() = default;

inline QPointF rotateAndNormalizePoint(double angle, const QPointF& vector)
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

std::optional<QPointF> segmentSegmentIntersection(QPointF A, QPointF B, QPointF C, QPointF D)
{
    if (intersect(A, B, C, D)
            || pointToLineSegmentDistanceSquared(A, C, D) == 0
            || pointToLineSegmentDistanceSquared(B, C, D) == 0
            || pointToLineSegmentDistanceSquared(C, A, B) == 0
            ||pointToLineSegmentDistanceSquared(D, A, B) == 0) {
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
    qSort(subtractedIntersectionPoints.begin(), subtractedIntersectionPoints.end(), QPointFLessThan);
    QVector<QPair<QPointF, int>> addedBackIntersectionPoints;
    for (const auto intersection: subtractedIntersectionPoints) {
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
    qDebug() << "-------------------------------------------------";
    qDebug() << "prunedInfill" << prunedInfill;
    qDebug() << "inset" << inset;
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

        qDebug() << "possiblePath" << possiblePath << polygonCircumferenceDistance(possiblePath);
        qDebug() << "reversePossiblePath" << possiblePathReverse << polygonCircumferenceDistance(possiblePathReverse);
        if (polygonCircumferenceDistance(possiblePath) <= polygonCircumferenceDistance(possiblePathReverse)) {
            if (connectedInfill.last() == possiblePath.first()) {
                possiblePath.pop_front();
            }
            connectedInfill << possiblePath;
            qDebug() << "first";
        } else {
            if (connectedInfill.last() == possiblePathReverse.first()) {
                possiblePathReverse.pop_front();
            }
            connectedInfill << possiblePathReverse;
            qDebug() << "second";
        }
        for (const auto pI: endConnectedInfill) {
            if (pI.first != connectedInfill.last()) {
                connectedInfill << pI.first;
            }
        }
    }
    return connectedInfill;
}

void InfillPattern::rectiLinearVerticalInfill()
{
    QVector<QPointF> infillPoints{};

    const auto outerBoundingBox = inputPolygon_.boundingRect();
    const auto startBoundary = outerBoundingBox.left();
    const auto endBoundary = outerBoundingBox.right();
    const auto top = outerBoundingBox.top();
    const auto bottom = outerBoundingBox.bottom();
    const auto insetTop = inset_.boundingRect().top();
    const auto insetBottom = inset_.boundingRect().bottom();
    const auto boundingPolygon = QPolygonF{{outerBoundingBox.topLeft(), outerBoundingBox.topRight(), outerBoundingBox.bottomRight(), outerBoundingBox.bottomLeft()}};

    auto up = true;

    for (auto x = startBoundary; x <= endBoundary; x += spacing_) {
        auto potentialTopPoint = QPointF{x, top};
        while((!boundingPolygon.containsPoint(potentialTopPoint, Qt::OddEvenFill) || qFuzzyCompare(potentialTopPoint.y(), insetTop)) && potentialTopPoint.y() <= bottom) {
            potentialTopPoint += QPointF{0, 1};
        }
        auto potentialBottomPoint = QPointF{x, bottom};
        while((!boundingPolygon.containsPoint(potentialBottomPoint, Qt::OddEvenFill) || qFuzzyCompare(potentialBottomPoint.y(), insetBottom)) && potentialBottomPoint.y() >= top) {
            potentialBottomPoint -= QPointF{0, 1};
        }

        if (potentialBottomPoint.y() - potentialTopPoint.y() > spacing_ &&
                boundingPolygon.containsPoint(potentialTopPoint, Qt::OddEvenFill) &&
                boundingPolygon.containsPoint(potentialBottomPoint, Qt::OddEvenFill)) {
            if (up) {
                infillPoints << potentialBottomPoint;
                infillPoints << potentialTopPoint;
            } else {
                infillPoints << potentialTopPoint;
                infillPoints << potentialBottomPoint;
            }
        }
        up = !up;
    }

    data_ << QVector<QVector<QPointF>>{infillPoints};
}

void InfillPattern::rectiLinearHorizontalInfill()
{
    QVector<QPointF> infillPoints{};

    const auto outerBoundingBox = inputPolygon_.boundingRect();
    const auto startBoundary = outerBoundingBox.top();
    const auto endBoundary = outerBoundingBox.bottom();
    const auto left = outerBoundingBox.left();
    const auto right = outerBoundingBox.right();
    auto directionRight = true;

    for (auto y = startBoundary; y <= endBoundary; y += spacing_) {
        auto potentialLeftPoint = QPointF{left, y};
        while(!inset_.containsPoint(potentialLeftPoint, Qt::OddEvenFill) && potentialLeftPoint.x() <= right) {
            potentialLeftPoint += QPointF{1, 0};
        }
        auto potentialRightPoint = QPointF{right, y};
        while(!inset_.containsPoint(potentialRightPoint, Qt::OddEvenFill) && potentialRightPoint.x() >= left) {
            potentialRightPoint -= QPointF{1, 0};
        }

        if (potentialRightPoint.x() - potentialLeftPoint.x() > spacing_ &&
                inset_.containsPoint(potentialLeftPoint, Qt::OddEvenFill) &&
                inset_.containsPoint(potentialRightPoint, Qt::OddEvenFill)) {
            if (directionRight) {
                infillPoints << potentialRightPoint;
                infillPoints << potentialLeftPoint;
            } else {
                infillPoints << potentialLeftPoint;
                infillPoints << potentialRightPoint;
            }
        }
        directionRight = !directionRight;
    }

    data_ << QVector<QVector<QPointF>>{infillPoints};
}

void InfillPattern::gridInfill()
{
    rectiLinearHorizontalInfill();
    const auto outerBoundingBox = inputPolygon_.boundingRect();
    const auto startBoundary = outerBoundingBox.left();
    const auto endBoundary = outerBoundingBox.right();
    const auto top = outerBoundingBox.top();
    const auto bottom = outerBoundingBox.bottom();
    const auto insetTop = inset_.boundingRect().top();
    const auto insetBottom = inset_.boundingRect().bottom();

    for (auto x = startBoundary; x <= endBoundary; x += spacing_) {
        auto potentialTopPoint = QPointF{x, top};
        while((!inset_.containsPoint(potentialTopPoint, Qt::OddEvenFill) || qFuzzyCompare(potentialTopPoint.y(), insetTop)) && potentialTopPoint.y() <= bottom) {
            potentialTopPoint += QPointF{0, 1};
        }
        auto potentialBottomPoint = QPointF{x, bottom};
        while((!inset_.containsPoint(potentialBottomPoint, Qt::OddEvenFill) || qFuzzyCompare(potentialBottomPoint.y(), insetBottom)) && potentialBottomPoint.y() >= top) {
            potentialBottomPoint -= QPointF{0, 1};
        }

        if (potentialBottomPoint.y() - potentialTopPoint.y() > spacing_ &&
                inset_.containsPoint(potentialTopPoint, Qt::OddEvenFill) &&
                inset_.containsPoint(potentialBottomPoint, Qt::OddEvenFill)) {
            QVector<QPointF> verticalLine;
            verticalLine << potentialBottomPoint << potentialTopPoint;
            data_ << verticalLine;
        }
    }
}

void InfillPattern::concentricInfill()
{
    const auto insetSpacing = spacing_ / 2;
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
        currentBorder = insetPolygon(currentBorder, insetSpacing);
    }
}

void InfillPattern::linearInfill()
{
    rectiLinearHorizontalInfill();
    const auto outerBoundingBox = inputPolygon_.boundingRect();
    const auto top = outerBoundingBox.top();
    const auto bottom = outerBoundingBox.bottom();
    const auto insetTop = inset_.boundingRect().top();
    const auto insetBottom = inset_.boundingRect().bottom();

    const auto startPositionX = outerBoundingBox.left() + outerBoundingBox.width() / 2;
    const auto startPositionY = outerBoundingBox.top() + outerBoundingBox.height() / 2;

    auto potentialTopPoint = QPointF{startPositionX, startPositionY};
    while((inset_.containsPoint(potentialTopPoint, Qt::OddEvenFill) || qFuzzyCompare(potentialTopPoint.y(), insetTop)) && potentialTopPoint.y() >= top) {
        potentialTopPoint += QPointF{0.7, -0.7};
    }
    potentialTopPoint -= QPointF{0.7, -0.7};
    auto potentialBottomPoint = QPointF{startPositionX, startPositionY};
    while((inset_.containsPoint(potentialBottomPoint, Qt::OddEvenFill)  || qFuzzyCompare(potentialBottomPoint.y(), insetBottom)) && potentialBottomPoint.y() <= bottom) {
        potentialBottomPoint += QPointF{-0.7, 0.7};
    }
    potentialBottomPoint -= QPointF{-0.7, 0.7};

    // TODO: Measure real length
    if (potentialBottomPoint.y() - potentialTopPoint.y() > spacing_ &&
            inset_.containsPoint(potentialTopPoint, Qt::OddEvenFill) &&
            inset_.containsPoint(potentialBottomPoint, Qt::OddEvenFill)) {
        QVector<QPointF> verticalLine;
        verticalLine << potentialBottomPoint << potentialTopPoint;
        data_ << verticalLine;
    }

    for (auto distance = spacing_; distance <= outerBoundingBox.width() / 2; distance += spacing_) {
        auto potentialRightTopPoint = QPointF{startPositionX + distance, startPositionY + distance * 0.7};
        while(inset_.containsPoint(potentialRightTopPoint, Qt::OddEvenFill) && potentialRightTopPoint.y() >= top) {
            potentialRightTopPoint += QPointF{0.7, -0.7};
        }
        potentialRightTopPoint -= QPointF{0.7, -0.7};
        auto potentialRightBottomPoint = QPointF{startPositionX + distance, startPositionY + distance * 0.7};
        while(inset_.containsPoint(potentialRightBottomPoint, Qt::OddEvenFill) && potentialRightBottomPoint.y() <= bottom) {
            potentialRightBottomPoint += QPointF{-0.7, 0.7};
        }
        potentialRightBottomPoint -= QPointF{-0.7, 0.7};

        // TODO: Measure real length
        if (potentialRightBottomPoint.y() - potentialRightTopPoint.y() > spacing_ &&
                inset_.containsPoint(potentialRightTopPoint, Qt::OddEvenFill) &&
                inset_.containsPoint(potentialRightBottomPoint, Qt::OddEvenFill)) {
            QVector<QPointF> verticalLine;
            verticalLine << potentialRightBottomPoint << potentialRightTopPoint;
            data_ << verticalLine;
        }

        auto potentialLeftTopPoint = QPointF{startPositionX - distance, startPositionY - distance * 0.7};
        while(inset_.containsPoint(potentialLeftTopPoint, Qt::OddEvenFill) && potentialLeftTopPoint.y() >= top) {
            potentialLeftTopPoint += QPointF{0.7, -0.7};
        }
        potentialLeftTopPoint -= QPointF{0.7, -0.7};
        auto potentialLeftBottomPoint = QPointF{startPositionX - distance, startPositionY - distance * 0.7};
        while(inset_.containsPoint(potentialLeftBottomPoint, Qt::OddEvenFill) && potentialLeftBottomPoint.y() <= bottom) {
            potentialLeftBottomPoint += QPointF{-0.7, 0.7};
        }
        potentialLeftBottomPoint -= QPointF{-0.7, 0.7};

        // TODO: Measure real length
        if (potentialLeftBottomPoint.y() - potentialLeftTopPoint.y() > spacing_ &&
                inset_.containsPoint(potentialLeftTopPoint, Qt::OddEvenFill) &&
                inset_.containsPoint(potentialLeftBottomPoint, Qt::OddEvenFill)) {
            QVector<QPointF> verticalLine;
            verticalLine << potentialLeftBottomPoint << potentialLeftTopPoint;
            data_ << verticalLine;
        }
    }

    for (auto distance = spacing_ / 2; distance <= outerBoundingBox.width() / 2; distance += spacing_) {
        auto potentialRightTopPoint = QPointF{startPositionX + distance, startPositionY - distance * 0.7};
        while(inset_.containsPoint(potentialRightTopPoint, Qt::OddEvenFill) && potentialRightTopPoint.y() >= top) {
            potentialRightTopPoint += QPointF{-0.7, -0.7};
        }
        potentialRightTopPoint -= QPointF{-0.7, -0.7};
        auto potentialRightBottomPoint = QPointF{startPositionX + distance, startPositionY - distance * 0.7};
        while(inset_.containsPoint(potentialRightBottomPoint, Qt::OddEvenFill) && potentialRightBottomPoint.y() <= bottom) {
            potentialRightBottomPoint += QPointF{0.7, 0.7};
        }
        potentialRightBottomPoint -= QPointF{0.7, 0.7};

        // TODO: Measure real length
        if (potentialRightBottomPoint.y() - potentialRightTopPoint.y() > spacing_ &&
                inset_.containsPoint(potentialRightTopPoint, Qt::OddEvenFill) &&
                inset_.containsPoint(potentialRightBottomPoint, Qt::OddEvenFill)) {
            QVector<QPointF> verticalLine;
            verticalLine << potentialRightBottomPoint << potentialRightTopPoint;
            data_ << verticalLine;
        }

        auto potentialLeftTopPoint = QPointF{startPositionX - distance, startPositionY + distance * 0.7};
        while(inset_.containsPoint(potentialLeftTopPoint, Qt::OddEvenFill) && potentialLeftTopPoint.y() >= top) {
            potentialLeftTopPoint += QPointF{-0.7, -0.7};
        }
        potentialLeftTopPoint -= QPointF{-0.7, -0.7};
        auto potentialLeftBottomPoint = QPointF{startPositionX - distance, startPositionY + distance * 0.7};
        while(inset_.containsPoint(potentialLeftBottomPoint, Qt::OddEvenFill) && potentialLeftBottomPoint.y() <= bottom) {
            potentialLeftBottomPoint += QPointF{0.7, 0.7};
        }
        potentialLeftBottomPoint -= QPointF{0.7, 0.7};

        // TODO: Measure real length
        if (potentialLeftBottomPoint.y() - potentialLeftTopPoint.y() > spacing_ &&
                inset_.containsPoint(potentialLeftTopPoint, Qt::OddEvenFill) &&
                inset_.containsPoint(potentialLeftBottomPoint, Qt::OddEvenFill)) {
            QVector<QPointF> verticalLine;
            verticalLine << potentialLeftBottomPoint << potentialLeftTopPoint;
            data_ << verticalLine;
        }
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
    const auto insetTop = inset_.boundingRect().top();
    const auto insetBottom = inset_.boundingRect().bottom();
    auto up = true;

    auto phaseIncrement = 1;
    for (auto x = startBoundary; x <= endBoundary; x += spacing_) {
        const auto decrementWith = std::min(0.1 * spacing_, 1.0);
        if (up) {
            for (auto y = bottom; y >= top; y -= decrementWith) {
                const auto potentialX = x + (spacing_ / 4) * sin((y - phaseIncrement * spacing_ / 4) * 4 / spacing_);
                auto potentialPoint = QPointF{potentialX, y};
                if (inset_.containsPoint(potentialPoint, Qt::OddEvenFill)  && !qFuzzyCompare(potentialPoint.y(), insetTop) && !qFuzzyCompare(potentialPoint.y(), insetBottom)) {
                    infillPoints << potentialPoint;
                }
            }
        } else {
            for (auto y = top; y <= bottom; y += decrementWith) {
                const auto potentialX = x + (spacing_ / 4) * sin((y - phaseIncrement * spacing_ / 4) * 4 / spacing_);
                auto potentialPoint = QPointF{potentialX, y};
                if (inset_.containsPoint(potentialPoint, Qt::OddEvenFill) && !qFuzzyCompare(potentialPoint.y(), insetTop) && !qFuzzyCompare(potentialPoint.y(), insetBottom)) {
                    infillPoints << potentialPoint;
                }
            }
        }

        up = !up;
        phaseIncrement = phaseIncrement + 1 % 2;
    }

    data_ << QVector<QVector<QPointF>>{infillPoints};
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
//            if (boundingPolygon.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
//            }
            currentPosition.setX(currentPosition.x() + L);
//            if (boundingPolygon.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
//            }
            currentPosition += rightDown;
//            if (boundingPolygon.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
//            }
            currentPosition.setX(currentPosition.x() + L);
//            if (boundingPolygon.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
//            }
            currentPosition += rightUp;
//            if (boundingPolygon.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
//            }
        }
        currentPosition = QPointF{currentPosition.x() + L, y + rightDown.y() * 2 + d};
        while (currentPosition.x() >= startBoundary) {
//            if (boundingPolygon.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
//            }
            currentPosition.setX(currentPosition.x() - L);
//            if (boundingPolygon.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
//            }
            currentPosition += leftUp;
//            if (boundingPolygon.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
//            }
            currentPosition.setX(currentPosition.x() - L);
//            if (boundingPolygon.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
//            }
            currentPosition += leftDown;
//            if (boundingPolygon.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
//            }
        }
        if (contiguousData.size() > 2) {
            contiguousData << contiguousData.first();
            infillPoints << contiguousData;
        }
    }

    data_ << infillPoints;
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
    return prunedInset;
}

void InfillPattern::pruneInfill(QPolygonF naiveInfill)
{
//    qDebug() << "-----------------------------------------------";
//    qDebug() << "naiveInfill" << naiveInfill;
//    QVector<QVector<QPair<QPointF, int>>> prunedInfill{};
//    auto outside = true;
//    for (int i = 0; i < naiveInfill.size() - 1; i++) {
//        QVector<QPair<QPointF, int>> connectedPrunedInfill{};
//        const auto firstPoint = naiveInfill[i];
//        const auto secondPoint = naiveInfill[i + 1];
//        const auto intersections = linePolygonIntersections(firstPoint, secondPoint, inset_);
//        if (intersections.isEmpty()) {
//            continue;
//        }
//        const auto sortedIntersections = sortPointsByDistance(firstPoint, intersections);
//        qDebug() << "intersections" << intersections;
//        qDebug() << "sortedIntersections" << sortedIntersections;
//        for (int j = 0; j < sortedIntersections.size() - 1; j += 2) {
////            if (outside) {
//                connectedPrunedInfill << sortedIntersections[j] << sortedIntersections[j + 1];
////            }
//            outside = !outside;
//        }
//        if (!connectedPrunedInfill.isEmpty()) {
//            prunedInfill << connectedPrunedInfill;
//            qDebug() << "added" << connectedPrunedInfill;
//        }
//    }
//    QVector<QVector<QPointF>> prunedInfill2{};
//    qDebug() << "apa1";
//    prunedInfill2 << QVector<QPointF>{};
//    if (!prunedInfill.empty()) {
//        QVector<QPointF> temp{};
//        for (const auto pI: prunedInfill[0]) {
//            temp << pI.first;
//        }
//        prunedInfill2[0] << temp;
//    }
//    qDebug()<< "prunedInfill" << prunedInfill;
//    for (int i = 0; i < prunedInfill.size() - 1; i++) {
//        const auto startConnectedInfill = prunedInfill[i];
//        const auto endConnectedInfill = prunedInfill[i + 1];
//        const auto startPointPair = startConnectedInfill.last();
//        const auto endPointPair = endConnectedInfill.first();
//        const auto startPoint = startPointPair.first;
//        const auto endPoint = endPointPair.first;
//        const auto insetStartIndex = startPointPair.second;
//        const auto insetEndIndex = endPointPair.second;
//        qDebug() << "start and end" << startPoint << endPoint << insetStartIndex << insetEndIndex;
//        if (false && insetStartIndex == insetEndIndex /*qFuzzyCompare(pointToLineSegmentDistanceSquared(endPoint, inset_[insetStartIndex], inset_[insetStartIndex + 1]), 0)*/) {
//            prunedInfill2[0] << startPoint << endPoint;
//        } else {

//            QPolygonF possiblePath = traverseFromStartToEnd(startPoint, endPoint, insetStartIndex, insetEndIndex, inset_);
//            const auto reversedPolygon = reversePolygonFOrientation(inset_, insetStartIndex, insetEndIndex);
//            qDebug() << "aaaaaaaaaaaaaaa" << inset_.size() << std::get<0>(reversedPolygon).size();
//            qDebug() << "aaaaaaaaaaaaaaa" << inset_ << std::get<0>(reversedPolygon);
//            const auto possiblePathReverse = traverseFromStartToEnd(startPoint, endPoint, std::abs(std::get<1>(reversedPolygon)), std::abs(std::get<2>(reversedPolygon)), std::get<0>(reversedPolygon));

//            qDebug() << "possiblePath" << possiblePath << polygonCircumferenceDistance(possiblePath);
//            qDebug() << "reversePossiblePath" << possiblePathReverse << polygonCircumferenceDistance(possiblePathReverse);
//            if (polygonCircumferenceDistance(possiblePath) <= polygonCircumferenceDistance(possiblePathReverse)) {
//                prunedInfill2[0] << possiblePath;
//                qDebug() << "first";
//            } else {
//                prunedInfill2[0] << possiblePathReverse;
//                qDebug() << "second";
//            }
//        }
////        for (int inset = insetStartIndex; inset < insetEndIndex - 1; inset++) {

////        }
//        QVector<QPointF> temp{};
//        for (const auto pI: endConnectedInfill) {
//            temp << pI.first;
//        }
//        prunedInfill2[0] << temp;
//    }
//    data_ << prunedInfill2;
}

void InfillPattern::update(QPolygonF &polygon)
{
    spacing_ = (100 / density_) * 12;
    data_.clear();
    inset_.clear();
    inset_ = insetPolygon(polygon, insetDistance_);
//    inset_ << inset_.first();
    inputPolygon_ = polygon;
    switch (pattern_) {
    case Pattern::None:
        break;
    case Pattern::Grid:
        gridInfill();
        break;
    case Pattern::RectiLinear: {
        rectiLinearVerticalInfill();
        const auto prunedInfill = ::pruneInfill(data_.first(), inset_);
        data_.clear();
        data_ << connectInfillAlongInset(prunedInfill, inset_, false);
            break;
        }
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

        const auto copy = data_;
        data_.clear();
        for (const auto& cc : copy) {
            const auto prunedInfill = ::pruneInfill(cc, inset_);
            data_ << connectInfillAlongInset(prunedInfill, inset_, true);
        }
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
    glColor3b(120, 100, 30);
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
