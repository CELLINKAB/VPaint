#include "InfillPattern.h"

#include <optional>

InfillPattern::InfillPattern() = default;

InfillPattern::~InfillPattern() = default;

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

    auto up = true;

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
    while((currentBorder.boundingRect().right() - currentBorder.boundingRect().left()) > insetSpacing &&
          (currentBorder.boundingRect().bottom() - currentBorder.boundingRect().top()) > insetSpacing) {
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

    for (auto y = top; y <= bottom; y+= (rightDown.y() + d) * 2) {
        QPolygonF contiguousData;
        QPointF currentPosition{startBoundary, y};
        while (currentPosition.x() <= endBoundary) {
            if (inset_.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
            }
            currentPosition.setX(currentPosition.x() + L);
            if (inset_.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
            }
            currentPosition += rightDown;
            if (inset_.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
            }
            currentPosition.setX(currentPosition.x() + L);
            if (inset_.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
            }
            currentPosition += rightUp;
            if (inset_.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
            }
        }
        currentPosition = QPointF{currentPosition.x() + L, y + rightDown.y() * 2 + d};
        while (currentPosition.x() >= startBoundary) {
            if (inset_.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
            }
            currentPosition.setX(currentPosition.x() - L);
            if (inset_.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
            }
            currentPosition += leftUp;
            if (inset_.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
            }
            currentPosition.setX(currentPosition.x() - L);
            if (inset_.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
            }
            currentPosition += leftDown;
            if (inset_.containsPoint(currentPosition, Qt::OddEvenFill)) {
                contiguousData << currentPosition;
            }
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

void InfillPattern::update(QPolygonF &polygon)
{
    spacing_ = (100 / density_) * 12;
    data_.clear();
    inset_.clear();
    inset_ = insetPolygon(polygon, insetDistance_);
    inputPolygon_ = polygon;
    switch (pattern_) {
    case Pattern::None:
        break;
    case Pattern::Grid:
        gridInfill();
        break;
    case Pattern::RectiLinear: {
        rectiLinearVerticalInfill();
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
        case Pattern::HoneyComb:
            honeycombInfill();
            break;
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
