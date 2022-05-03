#ifndef INFILLPATTERN_H
#define INFILLPATTERN_H

#include "VAC/vpaint_global.h"

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

#endif // INFILLPATTERN_H
