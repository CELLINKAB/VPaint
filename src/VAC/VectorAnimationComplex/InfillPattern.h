#ifndef INFILLPATTERN_H
#define INFILLPATTERN_H

#include "VAC/vpaint_global.h"

#include <optional>

// Creates and draws different infill patterns
// NOTE: Only works for convex shapes

class Q_VPAINT_EXPORT InfillPattern
{
public:
    enum class Pattern : uint8_t
    {
        None = 0,
        Grid = 1,
        HoneyComb = 2,
        RectiLinear = 3,
        Linear = 4,
        Concentric = 5,
        Gyroid = 6,
        MAX_VAL = 7
    };

    explicit InfillPattern();
    ~InfillPattern();

    void update(QPolygonF& polygon);

    [[nodiscard]]
    const QVector<QVector<QPointF>>& getPoints() const;

    void draw();

    [[nodiscard]] Pattern pattern() const;
    void setPattern(InfillPattern::Pattern pattern);
    [[nodiscard]] bool isApplyInfill() const;
    void setApplyInfill(bool apply);
    [[nodiscard]] int density() const;
    void setDensity(double density);

protected:
    QVector<QVector<QPointF>> data_{};
    Pattern pattern_{Pattern::None};
    QPolygonF inset_{};
    QVector<QPointF> normals_{};
    QPolygonF inputPolygon_;

private:
    void rectiLinearVerticalInfill();
    void rectiLinearHorizontalInfill();
    void gridInfill();
    void concentricInfill();
    void linearInfill();
    void gyroidInfill();
    void honeycombInfill();
    QPolygonF insetPolygon(QPolygonF polygon, double spacing);
    void pruneInfill(QPolygonF naiveInfill);
    bool applyInfill_{false};
    double density_{25};
    double spacing_ = 100;
    double insetDistance_ = 10;
};

Q_VPAINT_EXPORT std::optional<QPointF> segmentSegmentIntersection(QPointF A, QPointF B, QPointF C, QPointF D);
Q_VPAINT_EXPORT double pointToLineSegmentDistanceSquared(QPointF p, QPointF startPoint, QPointF endPoint);
Q_VPAINT_EXPORT double polygonCircumferenceDistance(QPolygonF polygon);
Q_VPAINT_EXPORT QPolygonF traverseFromStartToEnd(QPointF startPoint, QPointF endPoint, int insetStartIndex, int insetEndIndex, QPolygonF inset);
Q_VPAINT_EXPORT QVector<QPair<QPointF, int>> sortPointsByDistance(QPointF startPosition, QVector<QPair<QPointF, int>> foundIntersectionPoints);
Q_VPAINT_EXPORT QPolygonF reversePolygonFOrientation(QPolygonF polygon);
Q_VPAINT_EXPORT QVector<QVector<QPair<QPointF, int>>> pruneInfill(QPolygonF naiveInfill, QPolygonF inset);
Q_VPAINT_EXPORT QVector<QPointF> connectInfillAlongInset(QVector<QVector<QPair<QPointF, int>>> prunedInfill, QPolygonF inset);
Q_VPAINT_EXPORT std::tuple<int, int> startAndEndIndex(QPolygonF polygon, QPointF startPoint, QPointF endPoint);

#endif // INFILLPATTERN_H
