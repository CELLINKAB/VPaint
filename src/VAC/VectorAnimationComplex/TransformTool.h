// Copyright (C) 2012-2019 The VPaint Developers.
// See the COPYRIGHT file at the top-level directory of this distribution
// and at https://github.com/dalboris/vpaint/blob/master/COPYRIGHT
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef TRANSFORMTOOL_H
#define TRANSFORMTOOL_H

#include <QObject>

#include <QSet>
#include "../TimeDef.h"
#include "CellList.h"
#include "BoundingBox.h"

#include "Eigen.h"
#include "VAC/vpaint_global.h"

class ViewSettings;

namespace VectorAnimationComplex
{

class Q_VPAINT_EXPORT TransformTool: public QObject
{
    Q_OBJECT

public:
    // Constructor
    TransformTool(QObject * parent = 0);

    // Set state
    void setCells(const CellSet & cells);
    void setIdOffset(int idOffset);

    // Enum of all transform widgets IDs
    enum WidgetId
    {
        None,

        TopLeftScale,
        TopRightScale,
        BottomRightScale,
        BottomLeftScale,

        TopScale,
        RightScale,
        BottomScale,
        LeftScale,

        TopLeftRotate,
        TopRightRotate,
        BottomRightRotate,
        BottomLeftRotate,

        Pivot,

        MIN_WIDGET_ID = TopLeftScale,
        MAX_WIDGET_ID = Pivot
    };

    // Get which widget is currently hovered, if any
    WidgetId hovered() const;

    // Get pivot position
    Eigen::Vector2d pivotPosition(Time time) const;

    // Drawing
    void draw(const CellSet & cells, Time time, ViewSettings & viewSettings) const;

    // Picking
    void drawPick(const CellSet & cells, Time time, ViewSettings & viewSettings) const;
    void setHoveredObject(int id);
    void setNoHoveredObject();

    // Transform selection
    void beginTransform(double x0, double y0, Time time);
    // Added the angle parameter for the manual rotating
    bool continueTransform(double x, double y, double angle = 0.0);
    void endTransform();

    // Drag and drop transform tool
    void prepareDragAndDrop();
    void performDragAndDrop(double dx, double dy);
    void endDragAndDrop();

    bool setManualWidth(double newWidth, Time time);
    bool setManualHeight(double newHeight, Time time);
    bool setManualRotation(double angle, Time time);

private slots:
    void onKeyboardModifiersChanged();

private:
    // Disable copy and assignment
    TransformTool(const TransformTool &);
    TransformTool & operator=(const TransformTool &);

    CellSet cells_;
    int idOffset_;
    WidgetId hovered_;

    void glFillColor_(WidgetId id) const;
    void glStrokeColor_(WidgetId id) const;
    void glPickColor_(WidgetId id) const;

    void drawScaleWidget_(WidgetId id, const BoundingBox & bb, double size, ViewSettings & viewSettings) const;
    void drawPickScaleWidget_(WidgetId id, const BoundingBox & bb, double size, ViewSettings & viewSettings) const;

    void drawRotateWidget_(WidgetId id, const BoundingBox & bb, ViewSettings & viewSettings) const;
    void drawPickRotateWidget_(WidgetId id, const BoundingBox & bb, ViewSettings & viewSettings) const;

    void drawPivot_(const BoundingBox & bb, ViewSettings & viewSettings) const;
    void drawPickPivot_(const BoundingBox & bb, ViewSettings & viewSettings) const;

    // Pivot
    bool useAltTransform_() const;
    bool isPivotCached_  () const;

    Eigen::Vector2d manualPivotPosition_         () const;
    Eigen::Vector2d cachedTransformPivotPosition_() const;
    Eigen::Vector2d cachedPivotPosition_         () const;

    Eigen::Vector2d computePivotPosition_(Time time) const;
    Eigen::Vector2d computePivotPosition_(const BoundingBox & bb) const;

    Eigen::Vector2d pivotPosition_           (const BoundingBox & bb) const;
    Eigen::Vector2d noTransformPivotPosition_(const BoundingBox & bb) const;

    Eigen::Vector2d transformPivotPosition_       (WidgetId id, const BoundingBox & bb) const;
    Eigen::Vector2d defaultTransformPivotPosition_(WidgetId id, const BoundingBox & bb) const;
    Eigen::Vector2d altTransformPivotPosition_    (WidgetId id, const BoundingBox & bb) const;
    BoundingBox cellsGeometry(Time time) const;
    BoundingBox cellsGeometry() const;

    bool manualPivot_;
    bool draggingManualPivot_;
    double xManualPivot_, yManualPivot_;
    double xManualPivot0_, yManualPivot0_;

    bool dragAndDropping_;
    bool transforming_;
    bool rotating_;
    double xTransformPivot_, yTransformPivot_;
    double xTransformPivotAlt_, yTransformPivotAlt_;

    // Apply transformation
    bool isTransformConstrained_() const;
    KeyVertexSet draggedVertices_;
    KeyEdgeSet draggedEdges_;
    KeyFaceSet draggedFaces_;
    double x0_, y0_, dx_, dy_, x_, y_;
    BoundingBox bb0_, obb0_;
    double dTheta_;
    double dTheta0_;
    Time startTransformTime_;
};

}

#endif // TRANSFORMTOOL_H
