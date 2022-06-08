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

#include "Global.h"

#include "DevSettings.h"
#include "Scene.h"
#include "View.h"
#include "MainWindow.h"
#include "SettingsDialog.h"
#include "DevSettings.h"
#include "ColorSelector.h"
#include "Timeline.h"
#include "VectorAnimationComplex/KeyVertex.h"
#include "VAC/VectorAnimationComplex/EdgeGeometry.h"
#include "VAC/VectorAnimationComplex/EdgeSample.h"

#include <QDebug>
#include <QApplication>
#include <QKeyEvent>
#include <QHBoxLayout>
#include <QPushButton>
#include <QToolBar>
#include <QMenu>
#include <QAction>
#include <QSettings>
#include <QStatusBar>
#include <QDir>

namespace {
const constexpr auto SURFACE_ROUND_VERTICES = 100;
const constexpr auto IS_USED_STATUS_BAR = false;
}

// -------- Initialization --------

Global * global_ = 0;
Global * global() { return global_; }
void Global::initialize(MainWindow * w) { global_ = new Global(w); }

Global::Global(MainWindow * w) :
    toolMode_(SELECT),
    toolBar_(0),
    toolModeToolBar_(0),
    isScalingCorner_(false),
    isScalingEdge_(false),
    isRotating_(false),
    isDragAndDropping_(false),
    isDraggingPivot_(false),
    xSceneCursorPos_(0),
    ySceneCursorPos_(0),
    currentDisplayMode_(ILLUSTRATION),
    switchToDisplayMode_(OUTLINE),
    otherDisplayMode_(ILLUSTRATION_OUTLINE),
    mainWindow_(w),
    preferences_(),
    preferencesDialog_(0),
    settings_(0),
    documentDir_(QDir::home()),
    faceColor_(QColor::fromRgb(100, 100, 100, 10)),
    infillColor_(QColor::fromRgb(100, 100, 100, 10)),
    isDrawShapeFaceEnabled_(true),
    isShowAroundRectangleWhenDraw_(false),
    isShowVerticesOnSelection_(false),
    isShowGrid_(false),
    isGridSnapping_(true),
    highlightColorRatio_(1.2),
    highlightAlphaRatio_(2.0),
    selectColorRatio_(1.4),
    selectAlphaRatio_(3.0),
    surfaceScaleFactor_(1.0),
    gridSizeMM_(1.0),
    gridSize_(1.0),
    pasteDeltaX_(15),
    pasteDeltaY_(15),
    selectedRotation_(0.0),
    skippingCurveSamples_(2),
    selectedGeometry_{0.0, 0.0, 0.0, 0.0},
    mousePastePosition_{0.0, 0.0},
    sceneCenterPosition_{0.0, 0.0},
    arrowKeysMoveDelta_{1.0, 1.0},
    surfaceType_(VPaint::SurfaceType::None),
    surfaceBackGroundColor_{QColor(Qt::white)},
    surfaceBorderColor_{QColor(Qt::darkGray)},
    gridColor_{QColor(Qt::gray)},
    surfaceVertices_{},
    gridLines_{},
    gridHorizontalValues_{},
    gridVerticalValues_{}
{
    // Color selectors
    currentColor_ = new ColorSelector();
    currentColor_->setToolTip(tr("Current color (C)"));
    currentColor_->setStatusTip(tr("Click to open the color selector"));

    snapThreshold_ = new SpinBox();
    snapThreshold_->setCaption(tr(" snap threshold "));

    sculptRadius_ = new SpinBox();
    sculptRadius_->setCaption(tr(" sculpt radius "));

    // Event filter
    QCoreApplication::instance()->installEventFilter(this);

    // Status bar help
    statusBarHelp_ = new QLabel();
    statusBarHelp_->setText("Find help here.");
    w->statusBar()->addWidget(statusBarHelp_);
    connect(this, SIGNAL(keyboardModifiersChanged()), this, SLOT(updateStatusBarHelp()));
}

bool Global::deleteIsolatedVertices()
{
    return true;
}

bool Global::deleteShortEdges()
{
    return true;
}

Qt::KeyboardModifiers Global::keyboardModifiers()
{
    return keyboardModifiers_;
}

void Global::updateModifiers()
{
    Qt::KeyboardModifiers keyboardModifiers = QGuiApplication::queryKeyboardModifiers();
    if(keyboardModifiers_ != keyboardModifiers)
    {
        keyboardModifiers_ = keyboardModifiers;
        emit keyboardModifiersChanged();
    }
}

bool Global::eventFilter(QObject * /*watched*/, QEvent * event)
{
    // Every single event delivered by Qt go through this method first before
    // going to its target object, so keep it as lightweight as possible

    // It is used as a convenient way to fix a few event behaviours that were
    // not quite right out of the box.

    // --------------------- Detect modifier key presses --------------

    // Detect modifier key presses (Shift, Ctrl, Alt, etc.) and update application
    // state accordingly (e.g., indicate which modifiers are pressed in the status bar, or
    // redraw the scene, since highlighting color depends on which modifiers are pressed)

    // If a modifier is pressed or released, update the modifier state, and emit a signal
    // if this state has changed
    // If a modifier is pressed or released, update the modifier state, and emit a signal
    // if this state has changed
    if(event->type() == QEvent::KeyPress ||
       event->type() == QEvent::KeyRelease)
    {
        QKeyEvent * keyEvent = static_cast<QKeyEvent *>(event);
        if(keyEvent)
        {
            // Workaround for Mac delete key
            // This is needed because of a bug in QT 5 that has not been resolved as of 5.5.0
#ifdef Q_OS_MAC
            if(keyEvent->key() == Qt::Key_Backspace)
            {
                scene()->smartDelete();
            }
#endif
            if(keyEvent->key() == Qt::Key_Shift ||
               keyEvent->key() == Qt::Key_Alt ||
               keyEvent->key() == Qt::Key_Meta ||
               keyEvent->key() == Qt::Key_AltGr ||
               keyEvent->key() == Qt::Key_Control)
            {
                updateModifiers();
            }
        }

        // Continue normal processing of the event
        return false;
    }
    else if(event->type() == QEvent::FocusIn )
    {
        updateModifiers();

        // Continue normal processing of the event
        return false;
    }

    // --------------------- Resolve shortcut overloads --------------

    // Resolve shortcut overloads
    else if(event->type() == QEvent::Shortcut)
    {
        QShortcutEvent * shortcutEvent = static_cast<QShortcutEvent *>(event);

        if(shortcutEvent->isAmbiguous())
        {
            QKeySequence key = shortcutEvent->key();
            resolveAmbiguousShortcuts(key);

            // Stop processing of the event
            return true;

        }
        else
        {
            // Continue normal processing of the event
            return false;
        }
    }

    // --------------------- Keep standard behaviour --------------

    // Continue normal processing of the event
    return false;
}

void Global::resolveAmbiguousShortcuts(const QKeySequence & key)
{
    qDebug() << "Ambiguous shortcut:" << key;
}

Eigen::Vector2d Global::sceneCursorPos() const
{
    return Eigen::Vector2d(xSceneCursorPos_, ySceneCursorPos_);
}

void Global::setSceneCursorPos(const Eigen::Vector2d & pos)
{
    xSceneCursorPos_ = pos[0];
    ySceneCursorPos_ = pos[1];
}

QToolBar * Global::toolModeToolBar() const
{
    return toolModeToolBar_;
}

QToolBar * Global::toolBar() const
{
    return toolBar_;
}

void Global::createToolBars()
{
    // ----- Tool modes -----

    // Create toolbar
    toolBar_ = mainWindow()->addToolBar(tr("Toolbar"));
    toolBar_->setOrientation(Qt::Vertical);
    toolBar_->setMovable(false);
    mainWindow()->addToolBar(Qt::LeftToolBarArea, toolBar_);

    // Set toolbar size
    int iconWidth = 32;
    toolBar_->setIconSize(QSize(iconWidth,iconWidth));
    currentColor_->setIconSize(QSize(iconWidth,iconWidth));
    currentColor_->updateIcon();

    // Create actions (exclusive checkable)
    QActionGroup * actionGroup = new QActionGroup(this);
    for(int i=0; i<NUMBER_OF_TOOL_MODES; i++)
    {
        toolModeActions[i] = new ToolModeAction(static_cast<ToolMode>(i), actionGroup);
        toolModeActions[i]->setCheckable(true);
        toolModeActions[i]->setShortcutContext(Qt::ApplicationShortcut);
        toolBar_->addAction(toolModeActions[i]);
        connect(toolModeActions[i], SIGNAL(triggered(Global::ToolMode)),
                              this, SLOT(setToolMode(Global::ToolMode)));
    }

    // Select
    toolModeActions[SELECT]->setText(tr("Select and move (F1)"));
    toolModeActions[SELECT]->setIcon(QIcon(":/images/select.png"));
    toolModeActions[SELECT]->setStatusTip(tr("Select objects, move objects, glue objects together, and split curves."));
    toolModeActions[SELECT]->setShortcut(QKeySequence(Qt::Key_F1));

    // Sketch
    toolModeActions[SKETCH]->setText(tr("Sketch (F2)"));
    toolModeActions[SKETCH]->setIcon(QIcon(":/images/sketch.png"));
    toolModeActions[SKETCH]->setStatusTip(tr("Sketch curves."));
    toolModeActions[SKETCH]->setShortcut(QKeySequence(Qt::Key_F2));

    // Paint
    toolModeActions[PAINT]->setText(tr("Paint (F3)"));
    toolModeActions[PAINT]->setIcon(QIcon(":/images/paint.png"));
    toolModeActions[PAINT]->setStatusTip(tr("Paint an empty space or an existing object."));
    toolModeActions[PAINT]->setShortcut(QKeySequence(Qt::Key_F3));

    // Sculpt
    toolModeActions[SCULPT]->setText(tr("Sculpt (F4)"));
    toolModeActions[SCULPT]->setIcon(QIcon(":/images/sculpt.png"));
    toolModeActions[SCULPT]->setStatusTip(tr("Sculpt curves."));
    toolModeActions[SCULPT]->setShortcut(QKeySequence(Qt::Key_F4));

    // ----- Color selectors -----

    // Colors
    colorSelectorAction_ = toolBar_->addWidget(currentColor_);
    colorSelectorAction_->setText(tr("Color"));
    colorSelectorAction_->setToolTip(tr("Color (C)"));
    colorSelectorAction_->setStatusTip(tr("Click to open the color selector"));
    colorSelectorAction_->setShortcut(QKeySequence(Qt::Key_C));
    colorSelectorAction_->setShortcutContext(Qt::ApplicationShortcut);
    connect(colorSelectorAction_, SIGNAL(triggered()), currentColor_, SLOT(click()));

    // ----- Tool Options -----

    toolModeToolBar_ = new QToolBar("Action Bar");
    toolModeToolBar_->setIconSize(QSize(200,iconWidth));
    toolModeToolBar_->setMovable(false);
    mainWindow()->addToolBar(toolModeToolBar_);

    // ---------------------   Color   ------------------------

    actionChangeColor_ = new QAction(this);
    actionChangeColor_->setText(tr("Change color"));
    actionChangeColor_->setIcon(QIcon(":/images/change-color.png"));
    actionChangeColor_->setStatusTip(tr("Change the color of the selected cells"));
    //actionChangeColor_->setShortcut(QKeySequence(Qt::Key_C));
    actionChangeColor_->setShortcutContext(Qt::ApplicationShortcut);
    mainWindow()->addAction(actionChangeColor_);
    connect(actionChangeColor_, SIGNAL(triggered()), mainWindow()->scene(), SLOT(changeColor()));

    // ---------------------   Edges   ------------------------

    actionChangeEdgeWidth_ = new QAction(this);
    actionChangeEdgeWidth_->setText(tr("Change edge width (W)"));
    actionChangeEdgeWidth_->setIcon(QIcon(":/images/change-width.png"));
    actionChangeEdgeWidth_->setStatusTip(tr("Change the width of the selected edges"));
    actionChangeEdgeWidth_->setShortcut(QKeySequence(Qt::Key_W));
    actionChangeEdgeWidth_->setShortcutContext(Qt::ApplicationShortcut);
    mainWindow()->addAction(actionChangeEdgeWidth_);
    connect(actionChangeEdgeWidth_, SIGNAL(triggered()), mainWindow()->scene(), SLOT(changeEdgeWidth()));


    // ---------------------   Faces   ------------------------

    actionCreateFace_ = new QAction(this);
    actionCreateFace_->setText(tr("Create Face (F)"));
    actionCreateFace_->setIcon(QIcon(":/images/create-face.png"));
    actionCreateFace_->setStatusTip(tr("Create a face whose boundary is the selected edges"));
    actionCreateFace_->setShortcut(QKeySequence(Qt::Key_F));
    actionCreateFace_->setShortcutContext(Qt::ApplicationShortcut);
    mainWindow()->addAction(actionCreateFace_);
    connect(actionCreateFace_, SIGNAL(triggered()), mainWindow()->scene(), SLOT(createFace()));

    actionAddCycles_ = new QAction(this);
    actionAddCycles_->setText(tr("Add Holes (H)"));
    actionAddCycles_->setIcon(QIcon(":/images/add-cycles.png"));
    actionAddCycles_->setStatusTip(tr("Add holes to the selected face, whose boundaries are the selected edges"));
    actionAddCycles_->setShortcut(QKeySequence(Qt::Key_H));
    actionAddCycles_->setShortcutContext(Qt::ApplicationShortcut);
    mainWindow()->addAction(actionAddCycles_);
    connect(actionAddCycles_, SIGNAL(triggered()), mainWindow()->scene(), SLOT(addCyclesToFace()));

    actionRemoveCycles_ = new QAction(this);
    actionRemoveCycles_->setText(tr("Remove Holes (Ctrl+H)"));
    actionRemoveCycles_->setIcon(QIcon(":/images/remove-cycles.png"));
    actionRemoveCycles_->setStatusTip(tr("Remove holes from the selected face, whose boundaries are the selected edges"));
    actionRemoveCycles_->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_H));
    actionRemoveCycles_->setShortcutContext(Qt::ApplicationShortcut);
    mainWindow()->addAction(actionRemoveCycles_);
    connect(actionRemoveCycles_, SIGNAL(triggered()), mainWindow()->scene(), SLOT(removeCyclesFromFace()));

    // ---------------------   Topological operations   ------------------------

    actionGlue_ = new QAction(this);
    actionGlue_->setText(tr("Glue"));
    actionGlue_->setToolTip(tr("Glue (G)"));
    actionGlue_->setIcon(QIcon(":/images/glue.png"));
    actionGlue_->setStatusTip(tr("Glue two endpoints or two curves together"));
    actionGlue_->setShortcut(QKeySequence(Qt::Key_G));
    actionGlue_->setShortcutContext(Qt::ApplicationShortcut);
    mainWindow()->addAction(actionGlue_);
    connect(actionGlue_, SIGNAL(triggered()), mainWindow()->scene(), SLOT(glue()));

    actionUnglue_ = new QAction(this);
    actionUnglue_->setText(tr("Explode"));
    actionUnglue_->setToolTip(tr("Explode (E)"));
    actionUnglue_->setIcon(QIcon(":/images/unglue.png"));
    actionUnglue_->setStatusTip(tr("Duplicate the selected objects to disconnect adjacent curves and surfaces"));
    //actionUnglue_->setShortcut(QKeySequence(Qt::CTRL + Qt::Key_G));
    actionUnglue_->setShortcut(QKeySequence(Qt::Key_E));
    actionUnglue_->setShortcutContext(Qt::ApplicationShortcut);
    mainWindow()->addAction(actionUnglue_);
    connect(actionUnglue_, SIGNAL(triggered()), mainWindow()->scene(), SLOT(unglue()));

    actionUncut_ = new QAction(this);
    actionUncut_->setText(tr("Simplify"));
    actionUncut_->setToolTip(tr("Simplify (Backspace)"));
    actionUncut_->setIcon(QIcon(":/images/simplify.png"));
    actionUncut_->setStatusTip(tr("Simplify the selected objects, by merging curves and surfaces together"));
    actionUncut_->setShortcut(QKeySequence(Qt::Key_Backspace));
    actionUncut_->setShortcutContext(Qt::ApplicationShortcut);
    mainWindow()->addAction(actionUncut_);
    connect(actionUncut_, SIGNAL(triggered()), mainWindow()->scene(), SLOT(uncut()));

    // Desired icon size
    double sideLength = 40;

    // ---------------------   Shared actions/options   ------------------------

    // None

    // ---------------------   Select options   ------------------------

    toolModeToolBar_->addAction(actionChangeColor_);
    toolModeToolBar_->addAction(actionChangeEdgeWidth_);
    separatorSelect1_ = toolModeToolBar_->addSeparator();
    toolModeToolBar_->addAction(actionCreateFace_);
    toolModeToolBar_->addAction(actionAddCycles_);
    toolModeToolBar_->addAction(actionRemoveCycles_);
    separatorSelect2_ = toolModeToolBar_->addSeparator();

    toolModeToolBar_->widgetForAction(actionChangeColor_)->setFixedSize(sideLength,sideLength);
    toolModeToolBar_->widgetForAction(actionChangeEdgeWidth_)->setFixedSize(sideLength,sideLength);
    toolModeToolBar_->widgetForAction(actionCreateFace_)->setFixedSize(sideLength,sideLength);
    toolModeToolBar_->widgetForAction(actionAddCycles_)->setFixedSize(sideLength,sideLength);
    toolModeToolBar_->widgetForAction(actionRemoveCycles_)->setFixedSize(sideLength,sideLength);

    toolModeToolBar_->addAction(actionGlue_);
    toolModeToolBar_->addAction(actionUnglue_);
    toolModeToolBar_->addAction(actionUncut_);
    toolModeToolBar_->widgetForAction(actionGlue_)->setFixedSize(sideLength+20,sideLength);
    toolModeToolBar_->widgetForAction(actionUnglue_)->setFixedSize(sideLength+20,sideLength);
    toolModeToolBar_->widgetForAction(actionUncut_)->setFixedSize(sideLength+20,sideLength);


    // ---------------------   Sketch options   ------------------------

    // Tablet pressure
    actionUseTabletPressure_ = new QAction(this);
    actionUseTabletPressure_->setCheckable(true);
    actionUseTabletPressure_->setChecked(true);
    toolModeToolBar_->addAction(actionUseTabletPressure_);
    toolModeToolBar_->widgetForAction(actionUseTabletPressure_)->setFixedSize(sideLength,sideLength);
    actionUseTabletPressure_->setText(tr("Toggle stylus pressure"));
    actionUseTabletPressure_->setIcon(QIcon(":/images/pressure.png"));
    actionUseTabletPressure_->setStatusTip(tr("Enable or disable stylus pressure (only for users with a pen tablet)"));
    //actionUseTabletPressure_->setShortcut(QKeySequence(Qt::Key_Backspace));
    actionUseTabletPressure_->setShortcutContext(Qt::ApplicationShortcut);
    mainWindow()->addAction(actionUseTabletPressure_);
    connect(actionUseTabletPressure_, SIGNAL(triggered()), this, SLOT(toggleStylusPressure()));

    // Edge width
    edgeWidth_ = new SpinBox();
    edgeWidth_->setCaption(tr(" pen width "));
    edgeWidth_->setValue(settings().edgeWidth());
    actionEdgeWidth_ = toolModeToolBar_->addWidget(edgeWidth_);
    connect(edgeWidth_, SIGNAL(valueChanged(double)), this, SLOT(setEdgeWidth_(double)));

    // Separator
    separatorSketch1_ = toolModeToolBar_->addSeparator();

    // Planar map mode
    actionPlanarMapMode_ = new QAction(this);
    actionPlanarMapMode_->setCheckable(true);
    actionPlanarMapMode_->setChecked(false);
    toolModeToolBar_->addAction(actionPlanarMapMode_);
    toolModeToolBar_->widgetForAction(actionPlanarMapMode_)->setFixedSize(110,sideLength);
    actionPlanarMapMode_->setText(tr("Toggle intersections"));
    actionPlanarMapMode_->setIcon(QIcon(":/images/planar-map-on.png"));
    actionPlanarMapMode_->setStatusTip(tr("When intersections are enabled, the sketched curve automatically splits existing curves and surfaces."));
    //actionPlanarMapMode_->setShortcut(QKeySequence(Qt::Key_Backspace));
    actionPlanarMapMode_->setShortcutContext(Qt::ApplicationShortcut);
    mainWindow()->addAction(actionPlanarMapMode_);
    connect(actionPlanarMapMode_, SIGNAL(triggered()), this, SLOT(togglePlanarMapMode()));

    // Separator
    separatorSketch2_ = toolModeToolBar_->addSeparator();

    // Snapping
    actionSnapMode_ = new QAction(this);
    actionSnapMode_->setCheckable(true);
    actionSnapMode_->setChecked(false);
    toolModeToolBar_->addAction(actionSnapMode_);
    toolModeToolBar_->widgetForAction(actionSnapMode_)->setFixedSize(110,sideLength);
    actionSnapMode_->setText(tr("Toggle snapping"));
    actionSnapMode_->setIcon(QIcon(":/images/snapping-on.png"));
    actionSnapMode_->setStatusTip(tr("When snapping is enabled, the sketched curve is automatically glued to existing curves."));
    //actionSnapMode_->setShortcut(QKeySequence(Qt::Key_Backspace));
    actionSnapMode_->setShortcutContext(Qt::ApplicationShortcut);
    mainWindow()->addAction(actionSnapMode_);
    connect(actionSnapMode_, SIGNAL(triggered()), this, SLOT(toggleSnapping()));

    // Edge width
    actionSnapThreshold_ = toolModeToolBar_->addWidget(snapThreshold_);

    // ---------------------   Sculpt options   ------------------------

    actionSculptRadius_ = toolModeToolBar_->addWidget(sculptRadius_);

    // ---------------------   Cut options   ------------------------

    // Set default Tool Mode
    setToolMode(SKETCH);
}

Global::DisplayMode Global::displayMode() const
{
    return currentDisplayMode_;
}

View * Global::activeView() const
{
    return mainWindow()->activeView();
}

View * Global::hoveredView() const
{
    return mainWindow()->hoveredView();
}

Time Global::activeTime() const
{
    return activeView()->activeTime();
}

Timeline * Global::timeline() const
{
    return mainWindow()->timeline();
}

void Global::setDisplayMode(Global::DisplayMode mode)
{
    if(currentDisplayMode_ == mode)
    {
        // do nothing
    }
    else
    {
        currentDisplayMode_  = mode;
    }
}

bool Global::showCanvas() const
{
    return mainWindow()->isShowCanvasChecked();
}

void Global::togglePlanarMapMode()
{
    if(actionPlanarMapMode_->isChecked())
    {
        actionPlanarMapMode_->setText(tr("Disable intersections"));
        actionPlanarMapMode_->setIcon(QIcon(":/images/planar-map-on.png"));
    }
    else
    {
        actionPlanarMapMode_->setText(tr("Enable intersections"));
        actionPlanarMapMode_->setIcon(QIcon(":/images/planar-map-off.png"));
    }
}

void Global::toggleSnapping()
{
    if(actionSnapMode_->isChecked())
    {
        actionSnapMode_->setText(tr("Disable snapping"));
        actionSnapMode_->setIcon(QIcon(":/images/snapping-on.png"));
        actionSnapThreshold_->setEnabled(true);
    }
    else
    {
        actionSnapMode_->setText(tr("Enable snapping"));
        actionSnapMode_->setIcon(QIcon(":/images/snapping-off.png"));
        actionSnapThreshold_->setEnabled(false);
    }
}

void Global::toggleStylusPressure()
{
    // Nothing to do
}

void Global::setScalingCorner(bool b)
{
    isScalingCorner_ = b;
    updateStatusBarHelp();
}

void Global::setScalingEdge(bool b)
{
    isScalingEdge_ = b;
    updateStatusBarHelp();
}

void Global::setRotating(bool b)
{
    isRotating_ = b;
    updateStatusBarHelp();
}

void Global::setDragAndDropping(bool b)
{
    isDragAndDropping_ = b;
    updateStatusBarHelp();
}

void Global::setDraggingPivot(bool b)
{
    isDraggingPivot_ = b;
    updateStatusBarHelp();
}

Global::ToolMode Global::toolMode() const
{
    if(mainWindow()->isEditCanvasSizeVisible())
    {
        return EDIT_CANVAS_SIZE;
    }
    else
    {
        return toolMode_;
    }
}

void Global::setToolMode(Global::ToolMode mode)
{
    scene()->deselectAll();
    // Check consistency with action state
    if(!toolModeActions[mode]->isChecked())
        toolModeActions[mode]->setChecked(true);

    // Set member variable
    toolMode_ = mode;

    // Hide everything
    actionChangeColor_->setVisible(false);
    actionChangeEdgeWidth_->setVisible(false);
    actionCreateFace_->setVisible(false);
    actionAddCycles_->setVisible(false);
    actionRemoveCycles_->setVisible(false);

    toolModeToolBar_->removeAction(actionGlue_);
    toolModeToolBar_->removeAction(actionUnglue_);
    toolModeToolBar_->removeAction(actionUncut_);
    actionPlanarMapMode_->setVisible(false);
    actionSnapMode_->setVisible(false);
    actionUseTabletPressure_->setVisible(false);
    actionEdgeWidth_->setVisible(false);
    actionSnapThreshold_->setVisible(false);
    actionSculptRadius_->setVisible(false);
    separatorSelect1_->setVisible(false);
    separatorSelect2_->setVisible(false);
    separatorSketch1_->setVisible(false);
    separatorSketch2_->setVisible(false);

    // Desired icon size
    double sideLength = 40;

    // Show relevant
    switch(toolMode_)
    {
    case SELECT:
        actionChangeColor_->setVisible(true);
        actionChangeEdgeWidth_->setVisible(true);
        actionCreateFace_->setVisible(true);
        actionAddCycles_->setVisible(true);
        actionRemoveCycles_->setVisible(true);
        toolModeToolBar_->addAction(actionGlue_);
        toolModeToolBar_->addAction(actionUnglue_);
        toolModeToolBar_->addAction(actionUncut_);
        toolModeToolBar_->widgetForAction(actionGlue_)->setFixedSize(sideLength+20,sideLength);
        toolModeToolBar_->widgetForAction(actionUnglue_)->setFixedSize(sideLength+20,sideLength);
        toolModeToolBar_->widgetForAction(actionUncut_)->setFixedSize(sideLength+20,sideLength);
        separatorSelect1_->setVisible(true);
        separatorSelect2_->setVisible(true);
        break;
    case SKETCH:
        actionPlanarMapMode_->setVisible(true);
        actionSnapMode_->setVisible(true);
        actionUseTabletPressure_->setVisible(true);
        actionSnapThreshold_->setVisible(true);
        actionEdgeWidth_->setVisible(true);
        separatorSketch1_->setVisible(true);
        separatorSketch2_->setVisible(true);
        break;
    //Maybe in future we will need to change some properties for polygons
    case DRAW_LINE:
    case DRAW_RECTANGLE:
    case DRAW_CIRCLE:
    case DRAW_TRIANGLE:
    case DRAW_RHOMBUS:
    case DRAW_PENTAGON:
    case DRAW_HEXAGON:
    case DRAW_HEPTAGON:
    case DRAW_OCTAGON:
        actionPlanarMapMode_->setVisible(true);
        actionSnapMode_->setVisible(true);
        actionUseTabletPressure_->setVisible(true);
        actionSnapThreshold_->setVisible(true);
        actionEdgeWidth_->setVisible(true);
        separatorSketch1_->setVisible(true);
        separatorSketch2_->setVisible(true);
        break;
    case SCULPT:
        actionSculptRadius_->setVisible(true);
        break;
    case PAINT:
        break;
    default:
        break;
    }

    // Even when there's no icon, make it high enough
    toolModeToolBar_->setMinimumHeight(50);

    // Update help
    updateStatusBarHelp();

    // Update scene
    mainWindow()->update();
    mainWindow()->updatePicking();
}


void Global::updateStatusBarHelp()
{
    if (!IS_USED_STATUS_BAR)
        return;

    Qt::KeyboardModifiers keys = keyboardModifiers();
    bool isCtrlDown = keys & Qt::ControlModifier;
    bool isShiftDown = keys & Qt::ShiftModifier;
    bool isAltDown = keys & Qt::AltModifier;

    QString message;
    if(isCtrlDown || isShiftDown || isAltDown)
    {
        message += "[";
        if(isCtrlDown)
        {
            message += QString(ACTION_MODIFIER_NAME_SHORT).toUpper();
            if(isShiftDown || isAltDown)
                message += ",";
        }
        if(isShiftDown)
        {
            message += "SHIFT";
            if(isAltDown)
                message += ",";
        }

        if(isAltDown)
        {
            message += "ALT";
        }
        message += "] ";
    }

    if (isScalingCorner_)
    {
        if(!isShiftDown)
            message += "Hold SHIFT to preserve proportions. ";
        if(!isAltDown)
            message += "Hold ALT to scale relative to center/pivot. ";
    }
    else if (isScalingEdge_)
    {
        if(!isAltDown)
            message += "Hold ALT to scale relative to center/pivot. ";
    }
    else if (isRotating_)
    {
        if(!isShiftDown)
            message += "Hold SHIFT to rotate by 45° only. ";
        if(!isAltDown)
            message += "Hold ALT to rotate relative to opposite corner. ";
    }
    else if (isDragAndDropping_)
    {
        if(!isShiftDown)
            message += "Hold SHIFT to constrain translation along 45° axes. ";
    }
    else if (isDraggingPivot_)
    {
        if(!isShiftDown)
            message += "Hold SHIFT to snap to center and corners of bounding box. ";
    }
    else if(toolMode() == SELECT)
    {
        if(!isCtrlDown && !isShiftDown && !isAltDown) {
            message += "Click to select highlighted object. Click on background to deselect all. Hold " + QString(ACTION_MODIFIER_NAME_SHORT).toUpper() + ", SHIFT, or ALT for more actions.";
        }
        else if(isCtrlDown && !isShiftDown && !isAltDown) {
            message += "Click on curve to insert end point. Click on face to insert point-in-face.";
        }
        else if(!isCtrlDown && isShiftDown && !isAltDown) {
            message += "Click to add highlighted object to the selection. Hold also ALT for different action.";
        }
        else if(!isCtrlDown && !isShiftDown && isAltDown) {
            message += "Click to remove highlighted object from the selection. Hold also SHIFT for different action.";
        }
        else if(!isCtrlDown && isShiftDown && isAltDown) {
            message += "Click to select unselected objects, or deselect selected objects.";
        }
        else {
            message += "No action available for this combination of keyboard modifiers.";
        }
    }
    else if(toolMode() == SKETCH)
    {
        if(!isCtrlDown && !isShiftDown && !isAltDown) {
            message += "Hold left mouse button to draw a curve. " + QString(ACTION_MODIFIER_NAME_SHORT).toUpper() + ": Change pen width. ALT: Change snap threshold.";
        }
        else if(isCtrlDown && !isShiftDown && !isAltDown) {
            message += "Hold left mouse button to change pen width.";
        }
        else if(!isCtrlDown && !isShiftDown && isAltDown) {
            message += "Hold left mouse button to change snap threshold.";
        }
        else if(isCtrlDown && !isShiftDown && isAltDown) {
            message += "Hold left mouse button to change both pen width and snap threshold.";
        }
        else {
            message += "No action available for this combination of keyboard modifiers.";
        }
    }
    else if(toolMode() == PAINT)
    {
        if(!isCtrlDown && !isShiftDown && !isAltDown) {
            message += "Click on closed region delimited by curves to fill. Click on object to change color. Click on background to change background color.";
        }
        else {
            message += "No action available for this combination of keyboard modifiers.";
        }
    }
    else if(toolMode() == SCULPT)
    {
        if(!isCtrlDown && !isShiftDown && !isAltDown) {
            message += "Hold left mouse button (LMB) to drag endpoint, or drag curve within radius. " + QString(ACTION_MODIFIER_NAME_SHORT).toUpper() + ": radius. SHIFT: smooth. ALT: thickness.";
        }
        else if(isCtrlDown && !isShiftDown && !isAltDown) {
            message += "Hold LMB to change the radius of the sculpting tool. Note: radius not visible if cursor too far from curve.";
        }
        else if(!isCtrlDown && isShiftDown && !isAltDown) {
            message += "Hold LMB to smooth curve within radius.";
        }
        else if(!isCtrlDown && !isShiftDown && isAltDown) {
            message += "Hold LMB to change thickness of curve within radius. Trick: use a large radius to edit thickness of the whole curve.";
        }
        else {
            message += "No action available for this combination of keyboard modifiers.";
        }
    }

    statusBarHelp_->setText(message);
}

// Other getters
DevSettings * Global::devSettings() { return DevSettings::instance(); }
MainWindow * Global::mainWindow() const { return mainWindow_; }
Settings & Global::settings() { return preferences_; }
VPaint::Scene * Global::scene() const {return mainWindow()->scene();}

QColor Global::edgeColor()
{
    return currentColor_->color();
}

QColor Global::faceColor()
{
    return faceColor_;
}

QColor Global::infillColor()
{
    return infillColor_;
}

void Global::setEdgeColor(const QColor& newColor)
{
    if (newColor.isValid())
    {
        currentColor_->setColor(newColor);
        emit edgeColorChanged();
    }
}

void Global::setFaceColor(const QColor& newColor)
{
    if (newColor.isValid())
    {
        faceColor_ = newColor;
        emit faceColorChanged();
    }
}

void Global::setFaceAlpha(int alpha)
{
    faceColor_.setAlpha(alpha);
    emit faceColorChanged();
}

void Global::setInfillColor(const QColor &newColor)
{
    if (newColor.isValid())
    {
        infillColor_ = newColor;
        emit infillColorChanged();
    }
}

bool Global::isShowAroundRectangleWhenDraw() const
{
    return isShowAroundRectangleWhenDraw_;
}

void Global::setShowAroundRectangleWhenDraw(bool isShow)
{
    isShowAroundRectangleWhenDraw_ = isShow;
}

bool Global::isDrawShapeFaceEnabled() const
{
    return isDrawShapeFaceEnabled_;
}

void Global::setDrawShapeFaceEnabled(bool isEnabled)
{
    isDrawShapeFaceEnabled_ = isEnabled;
}

double Global::highlightColorRatio() const
{
    return highlightColorRatio_;
}

void Global::setHighlightColorRatio(double ratio)
{
    highlightColorRatio_ = ratio;
}

double Global::highlightAlphaRatio() const
{
    return highlightAlphaRatio_;
}

void Global::setHighlightAlphaRatio(double ratio)
{
    highlightAlphaRatio_ = ratio;
}

double Global::selectColorRatio() const
{
    return selectColorRatio_;
}

void Global::setSelectColorRatio(double ratio)
{
    selectColorRatio_ = ratio;
}

double Global::selectAlphaRatio() const
{
    return selectAlphaRatio_;
}

void Global::setSelectAlphaRatio(double ratio)
{
    selectAlphaRatio_ = ratio;
}

bool Global::isShowVerticesOnSelection() const
{
    return isShowVerticesOnSelection_;
}

void Global::setShowVerticesOnSelection(bool isShow)
{
    isShowVerticesOnSelection_ = isShow;
}

double Global::pasteDeltaX()
{
    return (isShowGrid_ && isGridSnapping_) ? gridSize_ : pasteDeltaX_;
}

double Global::pasteDeltaY()
{
    return (isShowGrid_ && isGridSnapping_) ? gridSize_ : pasteDeltaY_;
}

void Global::setPasteDelta(double dx, double dy)
{
    pasteDeltaX_ = dx;
    pasteDeltaY_ = dy;
}

void Global::setPasteDelta(double delta)
{
    pasteDeltaX_ = delta;
    pasteDeltaY_ = delta;
}

// Update the selection geometry for display the correct selection size/rotation in the Shape parameter bar
void Global::updateSelectedGeometry(double x, double y, double w, double h, double rotation)
{
    selectedRotation_ = rotation;
    updateSelectedGeometry(x, y, w, h);
}

void Global::updateSelectedGeometry(double x, double y, double w, double h)
{
    selectedGeometry_.setX(x);
    selectedGeometry_.setY(y);
    selectedGeometry_.setWidth(w);
    selectedGeometry_.setHeight(h);

    updateSelectedGeometry();
}

void Global::updateSelectedGeometry()
{
    emit interactiveGeometryChanged();
}

// Storing the mouse position on scene for be able to paste in right position
// Emiting rightMouseClicked() for show copy/paste popup
void Global::storeMousePastePos()
{
    mousePastePosition_.setX(xSceneCursorPos_);
    mousePastePosition_.setY(ySceneCursorPos_);
    emit rightMouseClicked();
}

void Global::setSurfaceType(VPaint::SurfaceType newSurfaceType)
{
    surfaceType_ = newSurfaceType;
    updateSurfaceVertices();
    updateGrid();
}

void Global::setShowGrid(bool showGrid)
{
    isShowGrid_ = showGrid;
    updateGrid();
    scene()->emitChanged();
}

void Global::setGridSnapping(bool gridSnapping)
{
    isGridSnapping_ = gridSnapping;
}

void Global::setSurfaceScaleFactor(double scaleFactor)
{
    surfaceScaleFactor_ = scaleFactor;
    gridSize_ = gridSizeMM_ * scaleFactor;
    sceneCenterPosition_.setX(scene()->left() + scene()->width() / 2);
    sceneCenterPosition_.setY(scene()->top() + scene()->height() / 2);
}

void Global::setGridSize(double gridSizeMM)
{
    gridSizeMM_ = gridSizeMM;
    gridSize_ = gridSizeMM * surfaceScaleFactor_;
    updateGrid();
    scene()->emitChanged();
}

void Global::setSurfaceColors(const QColor& bgColor, const QColor& borderColor, const QColor& gridColor)
{
    surfaceBackGroundColor_ = QColor(bgColor);
    surfaceBorderColor_ = QColor(borderColor);
    gridColor_ =  QColor(gridColor);
}

void Global::calculateSnappedHPosition(double& x)
{
    if (isShowGrid_ && isGridSnapping_)
    {
        auto xIndex = qCeil(round((scene()->left() + x) / gridSize_));
        xIndex = xIndex < gridHorizontalValues_.count() ? xIndex : gridHorizontalValues_.count() - 1;
        x = (gridHorizontalValues_.isEmpty() || xIndex < 0) ? x : gridHorizontalValues_[xIndex];
    }
}

void Global::calculateSnappedVPosition(double& y)
{
    if (isShowGrid_ && isGridSnapping_)
    {
        auto yIndex = qCeil(round((scene()->top() + y) / gridSize_));
        yIndex = yIndex < gridVerticalValues_.count() ? yIndex : gridVerticalValues_.count() - 1;
        y = (gridVerticalValues_.isEmpty() || yIndex < 0) ? y : gridVerticalValues_[yIndex];
    }
}

void Global::calculateSnappedPosition(double& x, double& y)
{
    if (isShowGrid_ && isGridSnapping_)
    {
        calculateSnappedHPosition(x);
        calculateSnappedVPosition(y);
    }
}

void Global::calculateSnappedPosition(QPointF& pos)
{
    if (isShowGrid_ && isGridSnapping_)
    {
        auto xRes = pos.x();
        auto yRes = pos.y();
        calculateSnappedPosition(xRes, yRes);
        pos.setX(xRes);
        pos.setY(yRes);
    }
}

QPointF Global::getSnappedPosition(double x, double y)
{
    auto xRes = x;
    auto yRes = y;
    calculateSnappedPosition(xRes, yRes);
    return QPointF(xRes, yRes);
}

QPointF Global::getSnappedPosition(const QPointF& pos)
{
    auto xRes = pos.x();
    auto yRes = pos.y();
    calculateSnappedPosition(xRes, yRes);
    return QPointF(xRes, yRes);
}

bool Global::isPointInSurface(const double x, const double y) const
{
    auto result = false;
    switch (surfaceType_) {
    case VPaint::SurfaceType::WellPlate:
    case VPaint::SurfaceType::PetriDish:
    {
        const auto radius = scene()->width() / 2;
        const auto deltaX = x - sceneCenterPosition_.x();
        const auto deltaY = y - sceneCenterPosition_.y();
        const auto distanceToCenter2 = deltaX * deltaX + deltaY * deltaY;
        result = distanceToCenter2 <= radius * radius;
        break;
    }
    case VPaint::SurfaceType::GlassSlide:
    {
        const auto xMin = scene()->left();
        const auto yMin = scene()->top();
        const auto xMax = xMin + scene()->width();
        const auto yMax = yMin + scene()->height();
        result = x >= xMin && y >= yMin && x <= xMax && y <= yMax;
        break;
    }
    default:
        break;
    }

    return result;
}

bool Global::isShapeInSurface(const VectorAnimationComplex::KeyVertexSet& vertices, const VectorAnimationComplex::KeyEdgeSet& edges, const QPointF delta) const
{
    // Check all vertices
    for (auto vertex : vertices)
    {
        const auto vx = vertex->pos()[0];
        const auto vy = vertex->pos()[1];
        if (!isPointInSurface(vx + delta.x(), vy + delta.y()))
            return false;
    }

    // Check samples of all edges(curves only)
    for (auto edge : edges)
    {
        if (edge->shapeType() == ShapeType::CURVE) {
            const auto samples = edge->geometry()->sampling();
            const auto inc = skippingCurveSamples() + 1;
            for (auto i = 0; i < samples.count(); i += inc)
            {
                const auto sample = samples[i];
                if (!isPointInSurface(sample[0] + delta.x(), sample[1] + delta.y()))
                    return false;
            }
        }
    }

    return true;
}

bool Global::isShapeInSurface(const QVector<QPointF>& vertices) const
{
    for (auto position : vertices)
    {
        if (!isPointInSurface(position.x(), position.y()))
            return false;
    }
    return true;
}

void Global::setArrowKeysMoveDelta(const double delta)
{
    arrowKeysMoveDelta_.setX(delta);
    arrowKeysMoveDelta_.setY(delta);
}

void Global::setArrowKeysMoveDelta(const double deltaX, const double deltaY)
{
    arrowKeysMoveDelta_.setX(deltaX);
    arrowKeysMoveDelta_.setY(deltaY);
}

void Global::setArrowKeysMoveDelta(const QPointF& delta)
{
    arrowKeysMoveDelta_.setX(delta.x());
    arrowKeysMoveDelta_.setY(delta.y());
}

void Global::setSkippingCurveSamples(const int value)
{
    skippingCurveSamples_ = value;
}

bool Global::useTabletPressure() const
{
    return actionUseTabletPressure_->isChecked();
}

double Global::edgeWidth() const
{
    return edgeWidth_->value();
}

void Global::setEdgeWidth_(double w)
{
    preferences_.setEdgeWidth(w);
}

void Global::updateSurfaceVertices()
{
    surfaceVertices_.clear();

    const auto x0 = scene()->left();
    const auto y0 = scene()->top();
    const auto w = scene()->width();
    const auto h = scene()->height();

    switch (surfaceType_) {
    case VPaint::SurfaceType::WellPlate:
    case VPaint::SurfaceType::PetriDish:
    {
        const auto x = x0 + w / 2;
        const auto y = y0 + h / 2;
        const auto r = w / 2;
        for (auto i = 0; i <= SURFACE_ROUND_VERTICES; i++)
        {
            auto angle = 2 * M_PI * i / SURFACE_ROUND_VERTICES;
            surfaceVertices_.append(QPointF(x0 + r * cos(angle) + x, y0 + r * sin(angle) + y));
        }
        break;
    }
    case VPaint::SurfaceType::GlassSlide:
    {
        surfaceVertices_.append(QPointF(x0, y0));
        surfaceVertices_.append(QPointF(x0 + w, y0));
        surfaceVertices_.append(QPointF(x0 + w, y0 + h));
        surfaceVertices_.append(QPointF(x0, y0 + h));
        break;
    }
    default:
        break;
    }
}

void Global::updateGrid()
{
    gridLines_.clear();
    gridHorizontalValues_.clear();
    gridVerticalValues_.clear();

    const auto x0 = scene()->left();
    const auto y0 = scene()->top();
    const auto w = scene()->width();
    const auto h = scene()->height();

    switch (surfaceType_) {
    case VPaint::SurfaceType::WellPlate:
    case VPaint::SurfaceType::PetriDish:
    {
        const auto x = x0 + w / 2;
        const auto y = y0 + h / 2;
        const auto rx = w / 2;
        const auto ry = h / 2;
        const auto gridLinesCount = qFloor(w / gridSize_);

        for (auto i = 0; i < gridLinesCount; i++)
        {
            auto lineX = x0 + i * gridSize_;
            auto lineY = y0 + i * gridSize_;
            auto deltaX = sqrt(abs(rx * rx - (lineY - y) * (lineY - y)));
            auto deltaY = sqrt(abs(ry * ry - (lineX - x) * (lineX - x)));

            gridLines_.append(QPair(QPointF(x - deltaX, lineY), QPointF(x + deltaX, lineY)));
            gridLines_.append(QPair(QPointF(lineX, y - deltaY), QPointF(lineX, y + deltaY)));
            gridHorizontalValues_.append(lineX);
            gridVerticalValues_.append(lineY);
        }
        break;
    }
    case VPaint::SurfaceType::GlassSlide:
    {
        const auto gridVLinesCount = qCeil(round(w / gridSize_));
        const auto gridHLinesCount = qCeil(round(h / gridSize_));

        for (auto i = 0; i < gridVLinesCount; i++)
        {
            auto lineX = x0 + i * gridSize_;
            gridLines_.append(QPair(QPointF(lineX, y0), QPointF(lineX, y0 + h)));
            gridHorizontalValues_.append(lineX);
        }

        for (auto i = 0; i < gridHLinesCount; i++)
        {
            auto lineY = y0 + i * gridSize_;
            gridLines_.append(QPair(QPointF(x0, lineY), QPointF(x0 + w, lineY)));
            gridVerticalValues_.append(lineY);
        }
        break;
    }
    default:
        break;
    }
}

void Global::setEdgeWidth(double w)
{
    if(edgeWidth() != w)
    {
        edgeWidth_->setValue(w);
    }

    preferences_.setEdgeWidth(w);
}

void Global::openPreferencesDialog()
{
    // Create preferences dialog
    if(!preferencesDialog_)
    {
        preferencesDialog_ = new SettingsDialog(mainWindow());
        connect(preferencesDialog_, SIGNAL(preferencesChanged()), this, SLOT(updateWidgetValuesFromPreferences()));
    }

    // Update and show references dialog
    preferencesDialog_->go();
}

void Global::updateWidgetValuesFromPreferences()
{
    edgeWidth_->setValue(preferences_.edgeWidth());
}

bool Global::planarMapMode() const
{
    return actionPlanarMapMode_->isChecked();
}

bool Global::snapMode() const
{
    return actionSnapMode_->isChecked();
}

double Global::snapThreshold() const
{
    return snapThreshold_->value();
}

void Global::setSnapThreshold(double newSnapThreshold)
{
    return snapThreshold_->setValue(newSnapThreshold);
}

double Global::sculptRadius() const
{
    return sculptRadius_->value();
}

void Global::setSculptRadius(double newRadius)
{
    return sculptRadius_->setValue(newRadius);
}

void Global::readSettings()
{
    QSettings qsettings;

    // Geometry of the window
    QSize size = qsettings.value("size", QSize(400, 400)).toSize();
    QPoint pos = qsettings.value("pos", QPoint(200, 200)).toPoint();
    mainWindow()->resize(size);
    mainWindow()->move(pos);

    // User settings
    settings().readFromDisk(qsettings);

    // Other settings
    snapThreshold_->setValue( qsettings.value("tools-sketch-snapthreshold", 15.0).toDouble() );
    sculptRadius_->setValue( qsettings.value("tools-sculpt-radius", 50.0).toDouble() );
}

void Global::writeSettings()
{
      QSettings qsettings;

      // Geometry of the window
      qsettings.setValue("size", mainWindow()->size());
      qsettings.setValue("pos", mainWindow()->pos());

      // User settings
      settings().writeToDisk(qsettings);

      // Other settings
      qsettings.setValue("tools-sketch-snapthreshold", snapThreshold_->value() );
      qsettings.setValue("tools-sculpt-radius", sculptRadius_->value() );
}

// ----- ToolModeAction -----

ToolModeAction::ToolModeAction(Global::ToolMode mode, QObject * parent) :
    QAction(parent),
    toolMode(mode)
{
    connect(this, SIGNAL(triggered()), this, SLOT(emitSpecializedTriggered()));
}

void ToolModeAction::emitSpecializedTriggered()
{
    emit triggered(toolMode);
}

void Global::setDocumentDir(const QDir & dir)
{
    documentDir_ = dir;
}

QDir Global::documentDir() const
{
    return documentDir_;
}
