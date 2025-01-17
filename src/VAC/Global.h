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

#ifndef GLOBAL_H
#define GLOBAL_H

// Define objects that we want accessible globally
//
// Example: global()->mainWindow()->update();
//          double w = global()->preferences().edgeWidth();

#include "Settings.h"

#include <QObject>
#include <QVector>
#include <QAction>
#include <QDoubleSpinBox>
#include "SpinBox.h"
#include "TimeDef.h"
#include <Eigen/Core>
#include <QLabel>
#include <QDir>
#include "VAC/vpaint_global.h"
#include "VectorAnimationComplex/CellList.h"

class QHBoxLayout;
class DevSettings;
class MainWindow;
class QToolBar;
class QAction;
class ToolModeAction;
class ColorSelector;
class QColor;
class SettingsDialog;
class View;
class Timeline;

namespace VPaint
{
class Scene;
enum class SurfaceType
{
    None = -1,
    PetriDish,
    WellPlate,
    GlassSlide
};
}

namespace VectorAnimationComplex
{
class VAC;
}

class Q_VPAINT_EXPORT Global: public QObject
{
    Q_OBJECT

public:
    // Initialization
    static void initialize(MainWindow * w);
    Global(MainWindow * w);

    // Tool Mode
    void createToolBars();
    enum ToolMode {
        // Used for array indexes, don't change the numbers!
        SELECT = 0,
        SKETCH,
        PAINT,
        SCULPT,
        //CUT,
        EDIT_CANVAS_SIZE, // This one is below "Number of tools" as it's not a mode interface-wise
        DRAW_LINE,
        DRAW_RECTANGLE,
        DRAW_CIRCLE,
        DRAW_TRIANGLE,
        DRAW_RHOMBUS,
        DRAW_PENTAGON,
        DRAW_HEXAGON,
        DRAW_HEPTAGON,
        DRAW_OCTAGON,
        NUMBER_OF_TOOL_MODES, // Keep this one last
    };
    ToolMode toolMode() const;

    // Menus
    void addSelectionActions(QMenu * selectionMenu);

    // Keyboard state
    Qt::KeyboardModifiers keyboardModifiers();

    // Tablet pressure
    bool useTabletPressure() const;

    // Edge width
    double edgeWidth() const;
    void setEdgeWidth(double w);

    // Planar map mode
    bool planarMapMode() const;

    // snapping
    bool snapMode() const;
    double snapThreshold() const;
    void setSnapThreshold(double newSnapThreshold);

    // Sculpting
    double sculptRadius() const;
    void setSculptRadius(double newRadius);

    // Automatic topological cleaning
    bool deleteIsolatedVertices();
    bool deleteShortEdges();

    // Cursor position
    Eigen::Vector2d sceneCursorPos() const;
    void setSceneCursorPos(const Eigen::Vector2d & pos);

    // Colors
    QColor edgeColor();
    QColor faceColor();
    QColor infillColor();

    void setEdgeColor(const QColor& newColor);
    void setFaceColor(const QColor& newColor);
    void setFaceAlpha(int alpha);
    void setInfillColor(const QColor& newColor);

    bool isShowAroundRectangleWhenDraw() const;
    void setShowAroundRectangleWhenDraw(bool isShow);

    bool isDrawShapeFaceEnabled() const;
    void setDrawShapeFaceEnabled(bool isEnabled);

    double highlightColorRatio() const;
    void setHighlightColorRatio(double ratio);

    double highlightAlphaRatio() const;
    void setHighlightAlphaRatio(double ratio);

    double selectColorRatio() const;
    void setSelectColorRatio(double ratio);

    double selectAlphaRatio() const;
    void setSelectAlphaRatio(double ratio);

    //For help and debug GCode generation
    bool isShowVerticesOnSelection() const;
    void setShowVerticesOnSelection(bool isShow);

    double pasteDeltaX();
    double pasteDeltaY();
    void setPasteDelta(double dx, double dy);
    void setPasteDelta(double delta);

    inline const QRectF& selectedGeometry() const { return selectedGeometry_; }
    inline double selectedRotation() const { return selectedRotation_; }
    void updateSelectedGeometry(double x, double y, double w, double h, double rotation);
    void updateSelectedGeometry(double x, double y, double w, double h);
    void updateSelectedGeometry();

    inline const QPointF& mousePastePosition() const { return mousePastePosition_; }
    void storeMousePastePos();

    inline VPaint::SurfaceType surfaceType() const { return surfaceType_; }
    void setSurfaceType(VPaint::SurfaceType newSurfaceType);

    inline bool isShowGrid() const { return isShowGrid_; }
    void setShowGrid(bool showGrid);

    inline double gridSizeMM() const { return gridSizeMM_; }
    inline double gridSize() const { return gridSize_; }
    void setGridSize(double gridSizeMM);

    inline bool isGridSnapping() const { return isShowGrid_ && isGridSnapping_; }
    void setGridSnapping(bool gridSnapping);

    inline double surfaceScaleFactor() const { return surfaceScaleFactor_; }
    void setSurfaceScaleFactor(double scaleFactor);

    inline const QColor& surfaceBackGroundColor() const { return surfaceBackGroundColor_; }
    inline const QColor& surfaceBorderColor() const { return surfaceBorderColor_; }
    inline const QColor& gridColor() const { return gridColor_; }

    void setSurfaceColors(const QColor& bgColor, const QColor& borderColor, const QColor& gridColor);

    inline const QList<QPointF>& surfaceVertices() const { return surfaceVertices_; }
    inline const QList<QPair<QPointF, QPointF>>& gridLines() const { return gridLines_; }

    void calculateSnappedHPosition(double& x);
    void calculateSnappedVPosition(double& y);
    void calculateSnappedPosition(double& x, double& y);
    void calculateSnappedPosition(QPointF& pos);
    QPointF getSnappedPosition(double x, double y);
    QPointF getSnappedPosition(const QPointF& pos);

    bool isPointInSurface(const double x, const double y) const;
    bool isShapeInSurface(const VectorAnimationComplex::KeyVertexSet& vertices, const VectorAnimationComplex::KeyEdgeSet& edges, const QPointF delta) const;
    bool isShapeInSurface(const VectorAnimationComplex::KeyVertexList& vertices, const VectorAnimationComplex::KeyEdgeList& edges) const;
    bool isShapeInSurface(const QVector<QPointF>& vertices) const;

    // A delta in millimeters used for moving shapes by
    // keyboard arrow keys(default: 1.0 mm for both directions)
    // If the grid snapping is enabled - delta will be the grid size
    inline const QPointF& arrowKeysMoveDelta() const { return arrowKeysMoveDelta_; }
    void setArrowKeysMoveDelta(const double delta);
    void setArrowKeysMoveDelta(const double deltaX, const double deltaY);
    void setArrowKeysMoveDelta(const QPointF& delta);

    // Count of skipping samples when checking they on inscribe into surface for
    // improve performance, because a curve line contains a lot of samples.
    inline int skippingCurveSamples() const { return skippingCurveSamples_; }
    void setSkippingCurveSamples(const int value);

    inline bool isDrawCircleAsCurve() const { return isDrawCircleAsCurve_; }
    void setIsDrawCircleAsCurve(const bool value);

    // Display modes
    enum DisplayMode {
        ILLUSTRATION,
        OUTLINE,
        ILLUSTRATION_OUTLINE
    };
    DisplayMode displayMode() const;
    void setDisplayMode(DisplayMode mode);
    bool showCanvas() const;

    // Active View and time
    View * activeView() const;
    View * hoveredView() const;
    Time activeTime() const;
    Timeline * timeline() const;

    // Other getters
    MainWindow * mainWindow() const;
    VPaint::Scene * scene() const;
    Settings & settings();
    DevSettings * devSettings();

    // Settings ( = user settings + application state )
    void readSettings();
    void writeSettings();

    // GUI elements owned by global
    QToolBar * toolModeToolBar() const;
    QToolBar * toolBar() const;

    // Directory from which paths in document are relative to
    void setDocumentDir(const QDir & dir);
    QDir documentDir() const;

signals:
    void keyboardModifiersChanged();
    void edgeColorChanged();
    void faceColorChanged();
    void infillColorChanged();
    void rightMouseClicked();
    void interactiveGeometryChanged();

public slots:
    void setToolMode(Global::ToolMode mode);
    void togglePlanarMapMode();
    void toggleSnapping();
    void toggleStylusPressure();

    void setScalingCorner(bool b);
    void setScalingEdge(bool b);
    void setRotating(bool b);
    void setDragAndDropping(bool b);
    void setDraggingPivot(bool b);

    // Open preference dialog
    void openPreferencesDialog();

    // Update widgets
    void updateWidgetValuesFromPreferences();

    // Help message
    void updateStatusBarHelp();

protected:
    // Global event filter
    bool eventFilter(QObject * watched, QEvent * event);

    // Update keyboard modifiers
    void updateModifiers();

    // Resolve ambiguous shortcuts
    void resolveAmbiguousShortcuts(const QKeySequence & key);

private slots:
    void setEdgeWidth_(double w);


private:
    // Tools
    ToolModeAction * toolModeActions [NUMBER_OF_TOOL_MODES];

    // Color selector
    QAction * colorSelectorAction_;

    // Tool Mode
    ToolMode toolMode_;
    QToolBar * toolBar_;

    // Tool options
    QToolBar * toolModeToolBar_;

    // Is a selection being transformed?
    bool isScalingCorner_;
    bool isScalingEdge_;
    bool isRotating_;
    bool isDragAndDropping_;
    bool isDraggingPivot_;

    // Select
    QAction * actionChangeColor_;
    QAction * actionChangeEdgeWidth_;
    QAction * actionCreateFace_;
    QAction * actionAddCycles_;
    QAction * actionRemoveCycles_;
    QAction * actionGlue_;
    QAction * actionUnglue_;
    QAction * actionUncut_;
    // Sketch
    QAction * actionPlanarMapMode_;
    QAction * actionSnapMode_;
    SpinBox * edgeWidth_;
    QAction * actionEdgeWidth_;
    SpinBox * snapThreshold_;
    QAction * actionSnapThreshold_;
    QAction * actionUseTabletPressure_;
    // Sculpt
    SpinBox * sculptRadius_;
    QAction * actionSculptRadius_;

    // Separators
    QAction * separatorSelect1_;
    QAction * separatorSelect2_;
    QAction * separatorSketch1_;
    QAction * separatorSketch2_;
    QAction * separatorSketch3_;

    // Scene cursor pos
    double xSceneCursorPos_, ySceneCursorPos_;

    // Colors
    ColorSelector * currentColor_;

    // Display modes
    DisplayMode currentDisplayMode_;
    DisplayMode switchToDisplayMode_;
    DisplayMode otherDisplayMode_;
    QAction * actionSwitchDisplayMode_;
    QAction * actionSwitchToOtherDisplayMode_;

    // Others
    MainWindow * mainWindow_;
    Settings preferences_;
    SettingsDialog * preferencesDialog_;
    DevSettings * settings_;
    Qt::KeyboardModifiers keyboardModifiers_;
    QDir documentDir_;

    // Status bar help
    QLabel * statusBarHelp_;

    // Colors
    QColor faceColor_;
    QColor infillColor_;

    bool isDrawShapeFaceEnabled_;
    bool isShowAroundRectangleWhenDraw_;
    bool isShowVerticesOnSelection_;
    bool isShowGrid_;
    bool isGridSnapping_;
    bool isDrawCircleAsCurve_;

    double highlightColorRatio_;
    double highlightAlphaRatio_;
    double selectColorRatio_;
    double selectAlphaRatio_;
    double surfaceScaleFactor_;
    double gridSizeMM_;
    double gridSize_;

    double pasteDeltaX_;
    double pasteDeltaY_;
    double selectedRotation_;
    double lastSurfaceHeight_;

    int skippingCurveSamples_;

    QRectF selectedGeometry_;
    QPointF mousePastePosition_;
    QPointF sceneCenterPosition_;

    QPointF arrowKeysMoveDelta_;

    VPaint::SurfaceType surfaceType_;
    QColor surfaceBackGroundColor_;
    QColor surfaceBorderColor_;
    QColor gridColor_;

    QList<QPointF> surfaceVertices_;

    QList<QPair<QPointF, QPointF>> gridLines_;
    QVector<double> gridHorizontalValues_;
    QVector<double> gridVerticalValues_;

    void updateSurfaceVertices();
    void updateGrid();
};

class Q_VPAINT_EXPORT ToolModeAction: public QAction
{
    Q_OBJECT

public:
    ToolModeAction(Global::ToolMode mode, QObject * parent = 0);

signals:
    void triggered(Global::ToolMode mode);

private slots:
    void emitSpecializedTriggered();

private:
    Global::ToolMode toolMode;
};

Q_VPAINT_EXPORT Global * global();

#endif // GLOBAL_H
