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

#include "View.h"

#include "Scene.h"
#include "Timeline.h"
#include "DevSettings.h"
#include "Global.h"
#include "Background/Background.h"
#include "Background/BackgroundRenderer.h"
#include "VectorAnimationComplex/VAC.h"
#include "VectorAnimationComplex/Cell.h"
#include "VectorAnimationComplex/KeyVertex.h"
#include "VectorAnimationComplex/CellList.h"
#include "VectorAnimationComplex/FaceCell.h"
#include "VectorAnimationComplex/KeyFace.h"

#include "Layer.h"
#include "MainWindow.h"

#include <QtDebug>
#include <QApplication>
#include <QPushButton>
#include <cmath>

// define mouse actions

#define  SELECT_ACTION                                      100
#define  ADDSELECT_ACTION                                   101
#define  DESELECT_ACTION                                    102
#define  TOGGLESELECT_ACTION                                103
#define  DESELECTALL_ACTION                                 104
#define  RECTANGLE_OF_SELECTION_ACTION                      105
#define  DRAG_AND_DROP_ACTION                               106
#define  SPLIT_ACTION                                       107
#define  TRANSFORM_SELECTION_ACTION                         108

#define  SKETCH_ACTION                                      200
#define  SKETCH_CHANGE_PEN_WIDTH_ACTION                     203
#define  SKETCH_CHANGE_SNAP_THRESHOLD_ACTION                204
#define  SKETCH_CHANGE_PEN_WIDTH_AND_SNAP_THRESHOLD_ACTION  205

#define  SCULPT_CHANGE_RADIUS_ACTION                        300
#define  SCULPT_DEFORM_ACTION                               301
#define  SCULPT_SMOOTH_ACTION                               302
#define  SCULPT_CHANGE_WIDTH_ACTION                         303

#define  PAINT_ACTION                                       400
#define  LINE_ACTION                                        401
#define  RECTANGLE_ACTION                                   402
#define  CIRCLE_ACTION                                      403
#define  TRIANGLE_ACTION                                    404
#define  RHOMBUS_ACTION                                     405
#define  PENTAGON_ACTION                                    406
#define  HEXAGON_ACTION                                     407
#define  HEPTAGON_ACTION                                    408
#define  OCTAGON_ACTION                                     409

namespace {
const constexpr auto CIRCLE_VERTICES = 40;
const constexpr auto POLYGON_ARROUND_VERTICES_FROM = 4;
const constexpr auto POLYGON_ARROUND_VERTICES_TO = 9;
const constexpr auto POLYGON_ARROUND_ALPHA = 0;
const constexpr auto POLYGON_ARROUND_LINE_SIZE = 0.5;
const constexpr auto DRAW_CIRCLES_AS_CURVES = false;
// Used to restrict the ability to draw very small shapes
const constexpr auto MIN_MOUSE_MOVE_TO_DRAW = 5;
}

View::View(VPaint::Scene * scene, QWidget * parent) :
    GLWidget(parent, true),
    scene_(scene),
    pickingImg_(0),
    pickingIsEnabled_(true),
    isMouseInSurface_(true),
    currentAction_(0),
    shapeStartX_(0),
    shapeStartY_(0),
    lastDragX_(0.0),
    lastDragY_(0.0),
    startDrawVertexCount(0),
    vac_(0)
{
    // View settings widget
    viewSettingsWidget_ = new ViewSettingsWidget(viewSettings_, this);
    connect(viewSettingsWidget_, SIGNAL(changed()), this, SLOT(update()));
    connect(viewSettingsWidget_, SIGNAL(changed()), this, SIGNAL(settingsChanged()));
    cameraTravellingIsEnabled_ = true;

    connect(this, SIGNAL(viewIsGoingToChange(int, int)), this, SLOT(updatePicking()));
    //connect(this, SIGNAL(viewIsGoingToChange(int, int)), this, SLOT(updateHighlightedObject(int, int)));
    connect(this, SIGNAL(viewIsGoingToChange(int, int)), this, SLOT(update()));

    connect(this, SIGNAL(viewIsBeingChanged(int, int)), this, SLOT(updateZoomFromView()));
    //connect(this, SIGNAL(viewIsBeingChanged(int, int)), this, SLOT(updatePicking()));
    //connect(this, SIGNAL(viewIsBeingChanged(int, int)), this, SLOT(updateHighlightedObject(int, int)));
    connect(this, SIGNAL(viewIsBeingChanged(int, int)), this, SLOT(update()));

    connect(this, SIGNAL(viewChanged(int, int)), this, SLOT(updateZoomFromView()));
    connect(this, SIGNAL(viewChanged(int, int)), this, SLOT(updatePicking()));
    connect(this, SIGNAL(viewChanged(int, int)), this, SLOT(updateHoveredObject(int, int)));
    connect(this, SIGNAL(viewChanged(int, int)), this, SLOT(update()));

    connect(global(), SIGNAL(keyboardModifiersChanged()), this, SLOT(handleNewKeyboardModifiers()));

    connect(global(), &Global::edgeColorChanged, this, [this]() { if(vac_) { vac_->changeEdgesColor(); }});
    connect(global(), &Global::faceColorChanged, this, [this]() { if(vac_) { vac_->changeFacesColor(); }});
    connect(global(), &Global::infillColorChanged, this, [this]() { if(vac_) { vac_->changeInfillColor(); }});
}

View::~View()
{
    deletePicking();
}

void View::initCamera()
{
    // Set 100% zoom and center canvas in view
    GLWidget_Camera2D camera;
    camera.setZoom(1.0);
    camera.setX(scene()->left() - 0.5*(scene()->width() - width()));
    camera.setY(scene()->top() - 0.5*(scene()->height() - height()));
    setCamera2D(camera);
}

VPaint::Scene * View::scene()
{
    return scene_;
}

void View::resizeEvent(QResizeEvent * event)
{
    if(autoCenterScene_)
        initCamera();

    GLWidget::resizeEvent(event);
}

void View::resizeGL(int width, int height)
{
    GLWidget::resizeGL(width, height);
    updatePicking();
}

void View::keyPressEvent(QKeyEvent *event)
{
    event->ignore();

    const auto moveDeltaX = global()->isGridSnapping() ?
                                global()->gridSizeMM() :
                                global()->arrowKeysMoveDelta().x();
    const auto moveDeltaY = global()->isGridSnapping() ?
                                global()->gridSizeMM() :
                                global()->arrowKeysMoveDelta().y();

    switch(event->key()){
    case Qt::Key_Up:
        moveSelected(0, -moveDeltaY);
        break;
    case Qt::Key_Down:
        moveSelected(0, moveDeltaY);
        break;
    case Qt::Key_Right:
        moveSelected(moveDeltaX, 0);
        break;
    case Qt::Key_Left:
        moveSelected(-moveDeltaX, 0);
        break;
    default:
        break;
    }
}


void View::keyReleaseEvent(QKeyEvent *event)
{
    event->ignore();
}

void View::handleNewKeyboardModifiers()
{
    vac_ = scene_->activeVAC();

    // Rectangle of selection
    if(vac_ && currentAction_ == RECTANGLE_OF_SELECTION_ACTION)
    {
        vac_->setSelectedCellsFromRectangleOfSelection(global()->keyboardModifiers());
    }

    // Update in any case, better be safe.
    emit allViewsNeedToUpdate();
}

MouseEvent View::mouseEvent() const
{
    MouseEvent me;
    me.x = mouse_Event_XScene_;
    me.y = mouse_Event_YScene_;
    me.left = mouse_LeftButton_;
    me.mid = mouse_MidButton_;
    me.right = mouse_RightButton_;
    me.alt = mouse_AltWasDown_;
    me.control = mouse_ControlWasDown_;
    me.shift = mouse_ShiftWasDown_;
    return me;
}

void View::update()
{
    GLWidget_Camera2D c = camera2D();
    c.setZoom(viewSettings_.zoom());
    setCamera2D(c);

    GLWidget::update();
}

void View::updateZoomFromView()
{
    viewSettings_.setZoom(zoom());
    viewSettingsWidget_->updateWidgetFromSettings();
    viewSettingsWidget_->updateSettingsFromWidgetSilent();
    GLWidget_Camera2D c = camera2D();
    c.setZoom(viewSettings_.zoom());
    setCamera2D(c);
}

int View::decideClicAction()
{
    vac_ = scene_->activeVAC();
    if (vac_)
    {
        if (global()->toolMode() == Global::SELECT && mouse_LeftButton_)
        {
            // Left = set selection
            if(!mouse_AltWasDown_ &&
               !mouse_ControlWasDown_ &&
               !mouse_ShiftWasDown_)
            {
                if (vac_->hoveredCell())
                    return SELECT_ACTION;
                else if (!vac_->hoveredTransformWidgetId())
                    return DESELECTALL_ACTION;
            }
            // Shift + Left = add to selection
            if(!mouse_AltWasDown_ &&
               !mouse_ControlWasDown_ &&
               mouse_ShiftWasDown_)
            {
                return ADDSELECT_ACTION;
            }
            // Alt + Left = remove from selection
            if(mouse_AltWasDown_ &&
               !mouse_ControlWasDown_ &&
               !mouse_ShiftWasDown_)
            {
                return DESELECT_ACTION;
            }
            // Alt + Shift + Left = toggle selection state
            if(mouse_AltWasDown_ &&
               !mouse_ControlWasDown_ &&
               mouse_ShiftWasDown_)
            {
                return TOGGLESELECT_ACTION;
            }
        }
    }

    return GLWidget::decideClicAction();
}

int View::decidePMRAction()
{
    vac_ = scene_->activeVAC();
    if (vac_)
    {
        // Selection
        if(global()->toolMode() == Global::SELECT)
        {
            // Left on cell
            if( vac_->hoveredCell() &&
                mouse_LeftButton_ &&
                !mouse_AltWasDown_ &&
                !mouse_ControlWasDown_ &&
                !mouse_ShiftWasDown_)
            {
                return DRAG_AND_DROP_ACTION;
            }
            // Left on transform widget
            else if( vac_->hoveredTransformWidgetId() &&
                     mouse_LeftButton_ &&
                     !mouse_ControlWasDown_)
            {
                return TRANSFORM_SELECTION_ACTION;
            }
            // Left on empty space
            else if (hoveredObject_.isNull() &&
                     mouse_LeftButton_ &&
                     !mouse_ControlWasDown_ )
            {
                return RECTANGLE_OF_SELECTION_ACTION;
            }
        }

        // Sketch
        else if(global()->toolMode() == Global::SKETCH)
        {
            // Left
            if(mouse_LeftButton_ &&
               !mouse_AltWasDown_ &&
               !mouse_ControlWasDown_ &&
               !mouse_ShiftWasDown_)
            {
                return SKETCH_ACTION;
            }
            // Ctrl + Left
            if(mouse_LeftButton_ &&
               !mouse_AltWasDown_ &&
               mouse_ControlWasDown_ &&
               !mouse_ShiftWasDown_)
            {
                return SKETCH_CHANGE_PEN_WIDTH_ACTION;
            }
            // Alt + Left
            if(mouse_LeftButton_ &&
               mouse_AltWasDown_ &&
               !mouse_ControlWasDown_ &&
               !mouse_ShiftWasDown_)
            {
                return SKETCH_CHANGE_SNAP_THRESHOLD_ACTION;
            }
            // Ctrl + Alt + Left
            if(mouse_LeftButton_ &&
               mouse_AltWasDown_ &&
               mouse_ControlWasDown_ &&
               !mouse_ShiftWasDown_)
            {
                return SKETCH_CHANGE_PEN_WIDTH_AND_SNAP_THRESHOLD_ACTION;
            }
        }

        // Sculpt
        else if(global()->toolMode() == Global::SCULPT)
        {
            // Left
            if(mouse_LeftButton_ &&
               !mouse_AltWasDown_ &&
               !mouse_ControlWasDown_ &&
               !mouse_ShiftWasDown_)
            {
                VectorAnimationComplex::Cell * hoveredCell =
                            vac_->hoveredCell();

                if(hoveredCell && hoveredCell->toVertexCell())
                {
                    return DRAG_AND_DROP_ACTION;
                }
                else
                {
                    return SCULPT_DEFORM_ACTION;
                }
            }
            // Ctrl + Left
            if(mouse_LeftButton_ &&
               !mouse_AltWasDown_ &&
               mouse_ControlWasDown_ &&
               !mouse_ShiftWasDown_)
            {
                return SCULPT_CHANGE_RADIUS_ACTION;
            }
            // Alt + Left
            if(mouse_LeftButton_ &&
               mouse_AltWasDown_ &&
               !mouse_ControlWasDown_ &&
               !mouse_ShiftWasDown_)
            {
                return SCULPT_CHANGE_WIDTH_ACTION;
            }
            // Shift + Left
            if(mouse_LeftButton_ &&
               !mouse_AltWasDown_ &&
               !mouse_ControlWasDown_ &&
               mouse_ShiftWasDown_)
            {
                return SCULPT_SMOOTH_ACTION;
            }
        }
        else if(global()->toolMode() == Global::DRAW_LINE)
        {
            return LINE_ACTION;
        }
        else if(global()->toolMode() == Global::DRAW_RECTANGLE)
        {
            return RECTANGLE_ACTION;
        }
            else if(global()->toolMode() == Global::DRAW_CIRCLE)
        {
            return CIRCLE_ACTION;
        }
        else if(global()->toolMode() == Global::DRAW_TRIANGLE)
        {
            return TRIANGLE_ACTION;
        }
        else if(global()->toolMode() == Global::DRAW_RHOMBUS)
        {
            return RHOMBUS_ACTION;
        }
        else if(global()->toolMode() == Global::DRAW_PENTAGON)
        {
            return PENTAGON_ACTION;
        }
        else if(global()->toolMode() == Global::DRAW_HEXAGON)
        {
            return HEXAGON_ACTION;
        }
        else if(global()->toolMode() == Global::DRAW_HEPTAGON)
        {
            return HEPTAGON_ACTION;
        }
        else if(global()->toolMode() == Global::DRAW_OCTAGON)
        {
            return OCTAGON_ACTION;
        }
    }

    return GLWidget::decidePMRAction();
}

int View::activeFrame() const
{
    return std::floor(viewSettings_.time().floatTime());
}

Time View::activeTime() const
{
    return viewSettings_.time();
}

void View::setActiveTime(Time t)
{
    viewSettings_.setTime(t);
    viewSettingsWidget_->updateWidgetFromSettings();
}

void View::setActive(bool isActive)
{
    viewSettingsWidget_->setActive(isActive);
}

void View::ClicEvent(int action, double x, double y)
{
    // It is View's responsibility to call update() or updatePicking()

    if(action==SPLIT_ACTION)
    {
        if(!hoveredObject_.isNull() || global()->toolMode() == Global::SKETCH)
        {
            vac_ = scene_->activeVAC();
            if(vac_)
            {
                vac_->split(x, y, interactiveTime(), true);

                emit allViewsNeedToUpdatePicking();
                updateHoveredObject(mouse_Event_X_, mouse_Event_Y_);
                emit allViewsNeedToUpdate();
            }
        }
    }
    else if(action==PAINT_ACTION)
    {
        Layer * layer = scene_->activeLayer();
        vac_ = layer ? layer->vac() : nullptr;
        if(vac_)
        {
            VectorAnimationComplex::Cell * paintedCell = vac_->paint(x, y, interactiveTime());
            if (!paintedCell)
            {
                layer->background()->setColor(global()->faceColor());
                scene_->emitChanged();
                scene_->emitCheckpoint();
            }

            emit allViewsNeedToUpdatePicking();
            updateHoveredObject(mouse_Event_X_, mouse_Event_Y_);
            emit allViewsNeedToUpdate();
        }
    }

    else if(action==SELECT_ACTION)
    {
        if(!hoveredObject_.isNull())
        {
            scene_->deselectAll();
            scene_->select(activeTime(),
                           hoveredObject_.index(),
                           hoveredObject_.id());
            emit allViewsNeedToUpdatePicking(); // required because selection bbox pickable
            updateHoveredObject(mouse_Event_X_, mouse_Event_Y_);
            emit allViewsNeedToUpdate();
        }
    }
    else if(action==DESELECTALL_ACTION)
    {
        scene_->deselectAll();
        emit allViewsNeedToUpdatePicking();
        updateHoveredObject(mouse_Event_X_, mouse_Event_Y_);
        emit allViewsNeedToUpdate();
    }
    else if(action==ADDSELECT_ACTION)
    {
        if(!hoveredObject_.isNull())
        {
            scene_->select(activeTime(),
                           hoveredObject_.index(),
                           hoveredObject_.id());

            emit allViewsNeedToUpdatePicking();
            updateHoveredObject(mouse_Event_X_, mouse_Event_Y_);
            emit allViewsNeedToUpdate();
        }
    }
    else if(action==DESELECT_ACTION)
    {
        if(!hoveredObject_.isNull())
        {
            scene_->deselect(activeTime(),
                             hoveredObject_.index(),
                             hoveredObject_.id());
            emit allViewsNeedToUpdatePicking();
            updateHoveredObject(mouse_Event_X_, mouse_Event_Y_);
            emit allViewsNeedToUpdate();
        }
    }
    else if(action==TOGGLESELECT_ACTION)
    {
        if(!hoveredObject_.isNull())
        {
            scene_->toggle(activeTime(),
                           hoveredObject_.index(),
                           hoveredObject_.id());
            emit allViewsNeedToUpdatePicking();
            updateHoveredObject(mouse_Event_X_, mouse_Event_Y_);
            emit allViewsNeedToUpdate();
        }
    }
    else
    {
        GLWidget::ClicEvent(action, x, y);
    }
}

void View::MoveEvent(double x, double y)
{
    // Boolean deciding if the scene must be redrawn even though only the mouse
    // has moved with no action performed. This is possible because depending on
    // where the mouse is, the action to-be-performed can be different, and
    // feedback to user on what is action would be must be given to user before
    // the action is undertaken
    bool mustRedraw = false;
    global()->setSceneCursorPos(Eigen::Vector2d(x,y));

    isMouseInSurface_ = global()->isPointInSurface(x, y);
    update();

    // Update highlighted object
    bool hoveredObjectChanged = updateHoveredObject(mouse_Event_X_, mouse_Event_Y_);
    if(hoveredObjectChanged)
        mustRedraw = true;

    // Update to-be-drawn straight line
    Qt::KeyboardModifiers keys = global()->keyboardModifiers();
    if( (global()->toolMode() == Global::SKETCH) )
    {
        if(keys & Qt::ControlModifier)
        {
            //scene_->updateCellsToConsiderForCutting();
            mustRedraw = true;
        }
        else
        {
            //scene_->resetCellsToConsiderForCutting();
            // Must be redrawn anyway to redraw the cursor
            mustRedraw = true;
        }
    }

    // Update to-be-sculpted edge
    if(global()->toolMode() == Global::SCULPT)
    {
        VectorAnimationComplex::VAC * vac = scene_->activeVAC();
        if (vac)
        {
            Time time = interactiveTime();
            vac->updateSculpt(x, y, time);
            mustRedraw = true;
        }
    }

    // Update to-be-painted face
    if(global()->toolMode() == Global::PAINT)
    {
        VectorAnimationComplex::VAC * vac = scene_->activeVAC();
        if (vac)
        {
            Time time = interactiveTime();
            vac->updateToBePaintedFace(x, y, time);
            mustRedraw = true;
        }
    }

    // Redraw if necessary
    if(mustRedraw)
    {
        // so that the highlighted object is also highlighted in other views
        // this is a matter of preference, we could call only "update()" if
        // we don't like this behaviour. But I like it, personally. Maybe
        // I could add it as a user preference
        emit allViewsNeedToUpdate();
    }
}

Time View::interactiveTime() const
{
    return viewSettings_.time();
}


void View::PMRPressEvent(int action, double x, double y)
{
    currentAction_ = action;

    // It is View's responsibility to call update() or updatePicking
    // for mouse PMR actions
    global()->setSceneCursorPos(Eigen::Vector2d(x,y));

    if (!isMouseInSurface_)
        return;

    if(action==SKETCH_ACTION)
    {
        // here, possibly,  the scene  has several layers  that it
        // knows about, as  well as which one is  active, and then
        // returns the active one.
        //
        // but the  scene doesn't know  at which time the  user is
        // drawing, since it depends  on the view.  (should active
        // layer depends on the view?  -> my current answer is no,
        // too  confusing.   But could  be  an option,  eventually
        // disable by default.  It would increases the flexibility
        // of the software).
        //
        // Current approach is then:
        //   1) The scene only knows which layer (ASG) is active
        //   2) The view only knows the time we are drawing in

        // Future ideas:
        //   Each  view  would  be  able  to see  the  scene  with
        //   different  translation/scale/rotation  (eg each  view
        //   has its own camera). Hence here, first the point
        //   (int xView, int yView) is converted into the point
        //   pos = (double xScene, double yScene)

        drawCurve(x, y, ShapeDrawPhase::DRAW_START);
    }
    else if(action==DRAG_AND_DROP_ACTION)
    {
        lastDragX_ = x;
        lastDragY_ = y;
        vac_->prepareDragAndDrop(mouse_PressEvent_XScene_, mouse_PressEvent_YScene_, interactiveTime());
    }
    else if(action==TRANSFORM_SELECTION_ACTION)
    {
        vac_->beginTransformSelection(mouse_PressEvent_XScene_, mouse_PressEvent_YScene_, interactiveTime());
    }
    else if(action==RECTANGLE_OF_SELECTION_ACTION)
    {
        vac_->beginRectangleOfSelection(x,y,interactiveTime());
    }
    else if(action==SCULPT_CHANGE_RADIUS_ACTION)
    {
        sculptStartRadius_ = global()->sculptRadius();
        sculptStartX_ = x;
        sculptStartY_ = y;
        sculptRadiusDx_ = 0;
        sculptRadiusDy_ = 0;

        //emit allViewsNeedToUpdatePicking();
        //updateHighlightedObject(mouse_Event_X_, mouse_Event_Y_);
        //emit allViewsNeedToUpdate();
    }
    else if(action==SKETCH_CHANGE_PEN_WIDTH_ACTION)
    {
        sculptStartRadius_ = global()->edgeWidth();
        sculptStartX_ = x;
        sculptStartY_ = y;
        sculptRadiusDx_ = 0;
        sculptRadiusDy_ = 0;

        //emit allViewsNeedToUpdatePicking();
        //updateHighlightedObject(mouse_Event_X_, mouse_Event_Y_);
        //emit allViewsNeedToUpdate();
    }
    else if(action==SKETCH_CHANGE_SNAP_THRESHOLD_ACTION)
    {
        sculptStartRadius_ = global()->snapThreshold();
        sculptStartX_ = x;
        sculptStartY_ = y;
        sculptRadiusDx_ = 0;
        sculptRadiusDy_ = 0;

        //emit allViewsNeedToUpdatePicking();
        //updateHighlightedObject(mouse_Event_X_, mouse_Event_Y_);
        //emit allViewsNeedToUpdate();
    }
    else if(action==SKETCH_CHANGE_PEN_WIDTH_AND_SNAP_THRESHOLD_ACTION)
    {
        sculptStartRadius_ = global()->edgeWidth();
        sculptStartRadius2_ = global()->snapThreshold();
        sculptStartX_ = x;
        sculptStartY_ = y;
        sculptRadiusDx_ = 0;
        sculptRadiusDy_ = 0;

        //emit allViewsNeedToUpdatePicking();
        //updateHighlightedObject(mouse_Event_X_, mouse_Event_Y_);
        //emit allViewsNeedToUpdate();
    }
    else if(action==SCULPT_DEFORM_ACTION)
    {
        sculptStartRadius_ = global()->sculptRadius();
        sculptStartX_ = x;
        sculptStartY_ = y;
        vac_->beginSculptDeform(x,y);

        //emit allViewsNeedToUpdatePicking();
        //updateHighlightedObject(mouse_Event_X_, mouse_Event_Y_);
        //emit allViewsNeedToUpdate();
    }
    else if(action==SCULPT_CHANGE_WIDTH_ACTION)
    {
        sculptStartRadius_ = global()->sculptRadius();
        sculptStartX_ = x;
        sculptStartY_ = y;
        vac_->beginSculptEdgeWidth(x,y);

        //emit allViewsNeedToUpdatePicking();
        //updateHighlightedObject(mouse_Event_X_, mouse_Event_Y_);
        //emit allViewsNeedToUpdate();
    }
    else if(action==SCULPT_SMOOTH_ACTION)
    {
        sculptStartRadius_ = global()->sculptRadius();
        sculptStartX_ = x;
        sculptStartY_ = y;
        vac_->beginSculptSmooth(x,y);

        //emit allViewsNeedToUpdatePicking();
        //updateHighlightedObject(mouse_Event_X_, mouse_Event_Y_);
        //emit allViewsNeedToUpdate();
    }
    else if(action == LINE_ACTION ||
            action == RECTANGLE_ACTION ||
            action == CIRCLE_ACTION ||
            action == TRIANGLE_ACTION ||
            action == RHOMBUS_ACTION ||
            action == PENTAGON_ACTION ||
            action == HEXAGON_ACTION ||
            action == HEPTAGON_ACTION ||
            action == OCTAGON_ACTION)
    {
        startDrawShape(x, y);
    }
    else
        GLWidget::PMRPressEvent(action, x, y);
}

void View::PMRMoveEvent(int action, double x, double y)
{
    global()->setSceneCursorPos(Eigen::Vector2d(x,y));

    if(action==DRAG_AND_DROP_ACTION)
    {
        // Check if a new position of the shape will be inside the surfaces area
        if (global()->isShapeInSurface(vac_->draggedVertices(), vac_->draggedEdges(), QPointF(x - lastDragX_, y - lastDragY_))) {
            vac_->performDragAndDrop(x, y);
            lastDragX_ = x;
            lastDragY_ = y;
        }
    }
    else if(action==TRANSFORM_SELECTION_ACTION)
    {
        vac_->continueTransformSelection(x, y);
    }
    else if(action == RECTANGLE_OF_SELECTION_ACTION)
    {
        vac_->continueRectangleOfSelection(x,y);
    }
    else if(action == SKETCH_ACTION)
    {
        drawCurve(x, y, ShapeDrawPhase::DRAW_PROCESS);
    }
    else if(action == LINE_ACTION)
    {
        drawLine(x, y, ShapeDrawPhase::DRAW_PROCESS);
    }
    else if(action == RECTANGLE_ACTION)
    {
        drawRectangle(x, y, ShapeDrawPhase::DRAW_PROCESS);
    }
    else if(action == CIRCLE_ACTION)
    {
        drawCircle(x, y, ShapeDrawPhase::DRAW_PROCESS);
    }
    else if(action == TRIANGLE_ACTION)
    {
        drawTriangle(x, y, ShapeDrawPhase::DRAW_PROCESS);
    }
    else if(action == RHOMBUS_ACTION)
    {
        drawPolygon(x, y, 4, 0, ShapeDrawPhase::DRAW_PROCESS);
    }
    else if(action == PENTAGON_ACTION)
    {
        drawPolygon(x, y, 5, -90, ShapeDrawPhase::DRAW_PROCESS);
    }
    else if(action == HEXAGON_ACTION)
    {
        drawPolygon(x, y, 6, 0, ShapeDrawPhase::DRAW_PROCESS);
    }
    else if(action == HEPTAGON_ACTION)
    {
        drawPolygon(x, y, 7, -90, ShapeDrawPhase::DRAW_PROCESS);
    }
    else if(action == OCTAGON_ACTION)
    {
        drawPolygon(x, y, 8, 22.5, ShapeDrawPhase::DRAW_PROCESS);
    }
    else if(action==SCULPT_CHANGE_RADIUS_ACTION)
    {
        // Increment how much we moved
        // method hiding the cursor and let it at the same position as on press
        // obviously can't work with pen tablets, since position is absolute...
        //sculptRadiusDx_ += x - sculptStartX_;
        //sculptRadiusDy_ += y - sculptStartY_;
        // hence just use the plain vanilla method
        sculptRadiusDx_ = x - sculptStartX_;
        sculptRadiusDy_ = y - sculptStartY_; // yes, this is useless, but can be useful later

        // Put mouse position back from where it was
        //QPoint p(mapToGlobal(QPoint(sculptStartX_,sculptStartY_)));
        //QCursor::setPos(p);

        // update radius
        double newRadius = sculptStartRadius_ + sculptRadiusDx_;
        if(newRadius < 0)
            newRadius *= -1;
        global()->setSculptRadius(newRadius);
    }
    else if(action==SKETCH_CHANGE_PEN_WIDTH_ACTION)
    {
        // Get delta
        sculptRadiusDx_ = x - sculptStartX_;
        sculptRadiusDy_ = y - sculptStartY_;

        // Get new radius
        double newRadius = sculptStartRadius_ + sculptRadiusDx_;
        if(newRadius < 0)
            newRadius *= -1;
        global()->setEdgeWidth(newRadius);

        // Prevent painted cursor gadget to move
        global()->setSceneCursorPos(Eigen::Vector2d(mouse_PressEvent_XScene_, mouse_PressEvent_YScene_));
    }
    else if(action==SKETCH_CHANGE_SNAP_THRESHOLD_ACTION)
    {
        // Get delta
        sculptRadiusDx_ = x - sculptStartX_;
        sculptRadiusDy_ = y - sculptStartY_;

        // Get new radius
        double newRadius = sculptStartRadius_ + sculptRadiusDx_;
        if(newRadius < 0)
            newRadius *= -1;
        global()->setSnapThreshold(newRadius);

        // Prevent painted cursor gadget to move
        global()->setSceneCursorPos(Eigen::Vector2d(mouse_PressEvent_XScene_, mouse_PressEvent_YScene_));
    }
    else if(action==SKETCH_CHANGE_PEN_WIDTH_AND_SNAP_THRESHOLD_ACTION)
    {
        // Get delta
        sculptRadiusDx_ = x - sculptStartX_;
        sculptRadiusDy_ = y - sculptStartY_;

        // Get new radius
        double newRadius = sculptStartRadius_ + sculptRadiusDx_;
        if(newRadius < 0)
            newRadius *= -1;
        global()->setEdgeWidth(newRadius);

        // Get new radius 2
        double newRadius2 = 0;
        if(sculptStartRadius_ > 0)
        {
            newRadius2 = sculptStartRadius2_ * newRadius / sculptStartRadius_;
        }
        else
        {
            newRadius2 = sculptStartRadius2_ + sculptRadiusDx_;
        }
        if(newRadius2 < 0)
            newRadius2 *= -1;
        global()->setSnapThreshold(newRadius2);

       // Prevent painted cursor gadget to move
        global()->setSceneCursorPos(Eigen::Vector2d(mouse_PressEvent_XScene_, mouse_PressEvent_YScene_));
    }
    else if(action==SCULPT_DEFORM_ACTION)
    {
        vac_->continueSculptDeform(x,y);
    }
    else if(action==SCULPT_CHANGE_WIDTH_ACTION)
    {
        vac_->continueSculptEdgeWidth(x,y);
    }
    else if(action==SCULPT_SMOOTH_ACTION)
    {
        vac_->continueSculptSmooth(x,y);
    }
    else
        GLWidget::PMRMoveEvent(action, x, y);

    emit allViewsNeedToUpdate();
    //updateView();
}

void View::PMRReleaseEvent(int action, double x, double y)
{
    currentAction_ = 0;

    global()->setSceneCursorPos(Eigen::Vector2d(x,y));

    if(action==DRAG_AND_DROP_ACTION)
    {
        vac_->completeDragAndDrop();
    }
    else if(action==TRANSFORM_SELECTION_ACTION)
    {
        vac_->endTransformSelection();
    }
    else if(action==RECTANGLE_OF_SELECTION_ACTION)
    {
        vac_->endRectangleOfSelection();
    }
    else if(action==SCULPT_CHANGE_RADIUS_ACTION)
    {
        vac_->updateSculpt(x, y, interactiveTime());
    }
    else if(action==SKETCH_CHANGE_PEN_WIDTH_AND_SNAP_THRESHOLD_ACTION)
    {

    }
    else if(action==SCULPT_DEFORM_ACTION)
    {
        vac_->endSculptDeform();
        vac_->updateSculpt(x, y, interactiveTime());
    }
    else if(action==SCULPT_CHANGE_WIDTH_ACTION)
    {
        vac_->endSculptEdgeWidth();
        vac_->updateSculpt(x, y, interactiveTime());
    }
    else if(action==SCULPT_SMOOTH_ACTION)
    {
        vac_->endSculptSmooth();
        vac_->updateSculpt(x, y, interactiveTime());
    }
    else if(action == SKETCH_ACTION)
    {
        drawCurve(x, y, ShapeDrawPhase::DRAW_END);
    }
    else if(action == LINE_ACTION)
    {
        drawLine(x, y, ShapeDrawPhase::DRAW_END);
    }
    else if(action == RECTANGLE_ACTION)
    {
        drawRectangle(x, y, ShapeDrawPhase::DRAW_END);
    }
    else if(action == CIRCLE_ACTION)
    {
        drawCircle(x, y, ShapeDrawPhase::DRAW_END);
    }
    else if(action == TRIANGLE_ACTION)
    {
        drawTriangle(x, y, ShapeDrawPhase::DRAW_END);
    }
    else if(action == RHOMBUS_ACTION)
    {
        drawPolygon(x, y, 4, 0, ShapeDrawPhase::DRAW_END);
    }
    else if(action == PENTAGON_ACTION)
    {
        drawPolygon(x, y, 5, -90, ShapeDrawPhase::DRAW_END);
    }
    else if(action == HEXAGON_ACTION)
    {
        drawPolygon(x, y, 6, 0, ShapeDrawPhase::DRAW_END);
    }
    else if(action == HEPTAGON_ACTION)
    {
        drawPolygon(x, y, 7, -90, ShapeDrawPhase::DRAW_END);
    }
    else if(action == OCTAGON_ACTION)
    {
        drawPolygon(x, y, 8, 22.5, ShapeDrawPhase::DRAW_END);
    }
    else
        GLWidget::PMRReleaseEvent(action, x, y);

    updateView();

}

/***********************************************************
 *              DRAWING
 */

void View::onBackgroundDestroyed_(Background * background)
{
    destroyBackgroundRenderer_(background);
}

BackgroundRenderer * View::getBackgroundRenderer_(Background * background)
{
    auto it = backgroundRenderers_.find(background);
    if (it != backgroundRenderers_.end())
    {
        return it.value();
    }
    else
    {
        return nullptr;
    }
}

BackgroundRenderer * View::createBackgroundRenderer_(Background * background)
{
    BackgroundRenderer * res = new BackgroundRenderer(background, this);
    connect(res, &BackgroundRenderer::backgroundDestroyed, this, &View::onBackgroundDestroyed_);
    backgroundRenderers_.insert(background, res);
    return res;
}

void View::destroyBackgroundRenderer_(Background * background)
{
    BackgroundRenderer * br = getBackgroundRenderer_(background);
    backgroundRenderers_.remove(background);
    delete br;
}

BackgroundRenderer * View::getOrCreateBackgroundRenderer_(Background * background)
{
    BackgroundRenderer * br = getBackgroundRenderer_(background);
    if (!br) {
        br = createBackgroundRenderer_(background);
    }
    return br;
}

void View::drawBackground_(Background* background, int frame)
{
    BackgroundRenderer* br = getOrCreateBackgroundRenderer_(background);

    switch (global()->surfaceType()) {
    case VPaint::SurfaceType::WellPlate:
    case VPaint::SurfaceType::PetriDish:
    {
        br->drawSurface(true);
        break;
    }
    case VPaint::SurfaceType::GlassSlide:
    {
        br->drawSurface();
        break;
    }
    default:
        br->draw(frame,
                 global()->showCanvas(),
                 scene_->left(), scene_->top(), scene_->width(), scene_->height(),
                 xSceneMin(), xSceneMax(), ySceneMin(), ySceneMax());
        break;
    }
}

void View::processRectangleOfSelection(double x, double y, ShapeDrawPhase drawPhase)
{
    if (!global()->isShowAroundRectangleWhenDraw())
        return;

    switch (drawPhase) {
    case ShapeDrawPhase::DRAW_START:
        vac_->beginRectangleOfSelection(x, y,interactiveTime());
        break;
    case ShapeDrawPhase::DRAW_PROCESS:
        vac_->continueRectangleOfSelection(x, y);
        break;
    case ShapeDrawPhase::DRAW_END:
        vac_->endRectangleOfSelection();
        break;
    default:
        break;
    }
}

void View::startDrawShape(double x, double y)
{
    processRectangleOfSelection(x, y, ShapeDrawPhase::DRAW_START);

    lastMousePos_ = QPoint(mouse_Event_X_, mouse_Event_Y_);

    QPointF pos = global()->getSnappedPosition(x, y);
    shapeStartX_ = pos.x();
    shapeStartY_ = pos.y();

    emit allViewsNeedToUpdate();
}

void View::endDrawShape()
{
    adjustCellsColors();    
    lastDrawnCells_.clear();
    vac_->endDrawShape();
    vac_->deselectAll();
    scene()->emitCheckpoint();
    emit scene()->changed();
}

void View::drawCurve(double x, double y, ShapeDrawPhase drawPhase)
{
    auto pos = global()->getSnappedPosition(x, y);
    auto xScene = pos.x();
    auto yScene = pos.y();

    double w = global()->settings().edgeWidth();
    if(mouse_isTablet_ &&  global()->useTabletPressure())
    {
        w *= 2 * mouse_tabletPressure_;
    }

    switch (drawPhase) {
    case ShapeDrawPhase::DRAW_START:
    {
        lastMousePos_ = QPoint(mouse_Event_X_, mouse_Event_Y_);
        vac_->beginSketchEdge(xScene, yScene, w, interactiveTime());
        emit allViewsNeedToUpdate();
        startDrawVertexCount = vac_->instantVertices().count();
        break;
    }
    case ShapeDrawPhase::DRAW_PROCESS:
    {
        const auto currentMousePos = QPoint(mouse_Event_X_, mouse_Event_Y_);
        if((currentMousePos - lastMousePos_).manhattanLength() < MIN_MOUSE_MOVE_TO_DRAW)
            return;

        if (global()->isPointInSurface(xScene, yScene)) {
            vac_->continueSketchEdge(xScene, yScene, w);
        } else {
            // Stop drawing a curve if mouse out of surface, because the curve will
            // be able to be visible outside the surface area, if continue moving the pressed
            // mouse outside the surface and return back to the surface area.
            forcedMouseRelease();
        }
        break;
    }
    case ShapeDrawPhase::DRAW_END:
    {
        vac_->endSketchEdge(false);
        if (vac_->instantVertices().count() - startDrawVertexCount < 2)
            break;

        lastDrawnCells_.clear();
        auto keyVertices = vac_->instantVertices();
        auto keyEdges = vac_->instantEdges();
        if (keyVertices.count() > 1 && !keyEdges.isEmpty()) {
            auto vertex1 = keyVertices.last();
            auto vertex2 = keyVertices.at(keyVertices.count() - 2);
            auto edge = keyEdges.last();
            vertex1->setShapeType(ShapeType::CURVE);
            vertex2->setShapeType(ShapeType::CURVE);
            edge->setShapeType(ShapeType::CURVE);
            lastDrawnCells_ << vertex1;
            lastDrawnCells_ << vertex2;
            lastDrawnCells_ <<  edge;
            if (global()->isDrawShapeFaceEnabled()) {
                for (auto cell : lastDrawnCells_)
                {
                    vac_->addToSelection(cell, false);
                }
                vac_->createFace(false);
                auto faceCell = vac_->faces().last();
                faceCell->setShapeType(ShapeType::CURVE);
                vac_->addToSelection(faceCell);
                lastDrawnCells_ << faceCell;
                // An one more edge has been created after the face has created
                // It connects the start and end vetices and we can ignore it
                auto lastEdge = vac_->instantEdges().last();
                lastEdge->setShapeType(ShapeType::CURVE);
                lastEdge->setIgnored(true);
                lastDrawnCells_ << lastEdge;
            }
            endDrawShape();
            scene()->emitShapeDrawn(ShapeType::CURVE);
        }
        break;
    }
    default:
        break;
    }
}

void View::drawLine(double x, double y, ShapeDrawPhase drawPhase)
{
    processRectangleOfSelection(x, y, drawPhase);

    switch (drawPhase) {
    case ShapeDrawPhase::DRAW_PROCESS:
    {
        drawShape(x, y, ShapeType::LINE, 2);
        break;
    }
    case ShapeDrawPhase::DRAW_END:
    {
        endDrawShape();
        scene()->emitShapeDrawn(ShapeType::LINE);
        break;
    }
    default:
        break;
    }
}

void View::drawCircle(double x, double y, ShapeDrawPhase drawPhase)
{
    processRectangleOfSelection(x, y, drawPhase);

    switch (drawPhase) {
    case ShapeDrawPhase::DRAW_PROCESS:
    {
        if (DRAW_CIRCLES_AS_CURVES) {
            // Draw circle using curved(will be enabled in future)
            // Now it's working fast but has some issues with infill
            drawShape(x, y, ShapeType::CIRCLE, CIRCLE_VERTICES);
        } else {
            // Draw circle as polygon, works fine but caused some
            // performance issues becasue it's used a lot of key vertices
            drawShape(x, y, ShapeType::POLYGON, CIRCLE_VERTICES, 0, true);
        }
        break;
    }
    case ShapeDrawPhase::DRAW_END:
    {
        endDrawShape();
        scene()->emitShapeDrawn(ShapeType::CIRCLE);
        break;
    }
    default:
        break;
    }
}

void View::drawTriangle(double x, double y, ShapeDrawPhase drawPhase)
{
    processRectangleOfSelection(x, y, drawPhase);

    switch (drawPhase) {
    case ShapeDrawPhase::DRAW_PROCESS:
    {
        drawShape(x, y, ShapeType::TRIANGLE, 3);
        break;
    }
    case ShapeDrawPhase::DRAW_END:
    {
        endDrawShape();
        scene()->emitShapeDrawn(ShapeType::TRIANGLE);
        break;
    }
    default:
        break;
    }
}

void View::drawRectangle(double x, double y, ShapeDrawPhase drawPhase)
{
    processRectangleOfSelection(x, y, drawPhase);

    switch (drawPhase) {
    case ShapeDrawPhase::DRAW_PROCESS:
    {
        drawShape(x, y, ShapeType::RECTANGLE, 4);
        break;
    }
    case ShapeDrawPhase::DRAW_END:
    {
         endDrawShape();
         scene()->emitShapeDrawn(ShapeType::RECTANGLE);
         break;
    }
    default:
        break;
    }
}

void View::drawPolygon(double x, double y, int countAngles, double rotation, ShapeDrawPhase drawPhase)
{
    processRectangleOfSelection(x, y, drawPhase);

    switch (drawPhase) {
    case ShapeDrawPhase::DRAW_PROCESS:
    {
        drawShape(x, y, ShapeType::POLYGON, countAngles, rotation);
        break;
    }
    case ShapeDrawPhase::DRAW_END:
    {
        // Since the polygon have to be inscribed in a circle/ellipse, so a transparent
        // bounding circle/ellipse is drawn around the polygons with more than 4 vertices.
        // This is used to correct the behavior of polygons, namely:
        // - displaying the same polygon size at the end of the drawing and after selection on the shape parameters bar;
        // - ability to draw a regular(equilateral) polygon(when holding down the "Shift" key during draw)
        // - correcting behavior of the grid snapping feature for polygons(works like as for circle/ellipse).
        // All cells of this bounding circle/ellipse are marked as ignored and aren't used for generate a G-code
        // It is now visible for debugging/demonstration, later it will be setted to transparent.
        if (countAngles > POLYGON_ARROUND_VERTICES_FROM && countAngles < POLYGON_ARROUND_VERTICES_TO)
        {
            const auto angleStep = M_PI / 50;
            QColor boundingColor = global()->edgeColor();
            boundingColor.setAlpha(POLYGON_ARROUND_ALPHA);

            const auto rx = global()->selectedGeometry().width() / 2;
            const auto ry = global()->selectedGeometry().height() / 2;
            const auto xCenter = global()->selectedGeometry().x() + rx;
            const auto yCenter = global()->selectedGeometry().y() + ry;
            vac_->beginSketchEdge(xCenter + rx, yCenter, POLYGON_ARROUND_LINE_SIZE, interactiveTime());
            for (auto angle = 0.0; angle <= M_PI * 2; angle += angleStep)
            {
                const auto newX = xCenter + rx * cos(angle);
                const auto newY = yCenter + ry * sin(angle);
                vac_->continueSketchEdge(newX, newY, POLYGON_ARROUND_LINE_SIZE);
            }
            vac_->endSketchEdge(false);

            auto keyEdge = vac_->instantEdges().last();
            keyEdge->setIgnored(true);
            vac_->assignShapeID(keyEdge);
            keyEdge->setColor(boundingColor);
            keyEdge->setHighlightedColor(boundingColor);
            keyEdge->setSelectedColor(boundingColor);
        }
        endDrawShape();
        scene()->emitShapeDrawn(ShapeType::POLYGON);
        break;
    }
    default:
        break;
    }
}

void View::adjustCellsColors()
{
    for (auto cell : lastDrawnCells_)
    {
        vac_->assignShapeID(cell);
        vac_->adjustSelectColors(cell);
    }
}

void View::drawShape(double x, double y, ShapeType shapeType, int countAngles, double initialRotation, bool drawingCircle)
{
    auto isSuccess = false;
    const auto currentMousePos = QPoint(mouse_Event_X_, mouse_Event_Y_);

    const auto pos = global()->getSnappedPosition(x, y);
    const auto xScene = pos.x();
    const auto yScene = pos.y();

    if((currentMousePos - lastMousePos_).manhattanLength() < MIN_MOUSE_MOVE_TO_DRAW)
        return;

    using CellSet = VectorAnimationComplex::CellSet;
    using Vertex = VectorAnimationComplex::KeyVertex;
    using Edge = VectorAnimationComplex::KeyEdge;
    using EggeSet = VectorAnimationComplex::KeyEdgeSet;
    using Cycle = VectorAnimationComplex::Cycle;

    QVector<QPointF> verticesPoints(countAngles);

    auto w = global()->settings().edgeWidth();
    if(mouse_isTablet_ &&  global()->useTabletPressure())
    {
        w *= 2 * mouse_tabletPressure_; // 2 so that a half-pressure would get the default width
    }

    auto leftX = xScene > shapeStartX_ ? shapeStartX_ : xScene;
    auto rightX = leftX == shapeStartX_ ? xScene : shapeStartX_;
    auto topY = yScene > shapeStartY_ ? shapeStartY_ : yScene;
    auto bottomY = topY == shapeStartY_ ? yScene : shapeStartY_;

    auto shapeWidth = rightX - leftX;
    auto shapeHeight = bottomY - topY;

    if (shapeType != ShapeType::LINE && global()->keyboardModifiers() == Qt::ShiftModifier)
    {
        if (shapeWidth > shapeHeight)
        {
            bottomY = topY + shapeWidth;
            shapeHeight = shapeWidth;
        }
        else if (shapeHeight > shapeWidth)
        {
            rightX = leftX + shapeHeight;
            shapeWidth = shapeHeight;
        }
    }

    auto clearLastCells = [this]()
    {
        if (!lastDrawnCells_.isEmpty())
        {
            vac_->deleteCells(lastDrawnCells_);
            lastDrawnCells_.clear();
        }
    };

    auto processDrawShape = [this, &verticesPoints, &w, shapeType, &drawingCircle, &clearLastCells, &isSuccess](bool isClosedShape = true)
    {
        const auto verticesCount = verticesPoints.count();

        // Check if a new shape is inside the surfaces area
        if (!global()->isShapeInSurface(verticesPoints))
            return;

        clearLastCells();

        const auto actualShapeType = drawingCircle ? ShapeType::CIRCLE : shapeType;

        //Draw Vertices
        QVector<Vertex*> vertices;
        for (auto point : verticesPoints)
        {
            auto vertex = vac_->newKeyVertex(interactiveTime(), Eigen::Vector2d(point.x(), point.y()));
            vertex->setShapeType(actualShapeType);
            vertex->setColor(global()->edgeColor());

            vertices << vertex;
            lastDrawnCells_ << vertex;
        }

        //Draw Edges
        EggeSet edges;
        for (auto i = 0; i < verticesCount - 1; i++)
        {
            auto keyEdge = vac_->newKeyEdge(interactiveTime(), vertices[i], vertices[i + 1], 0, w);
            keyEdge->setShapeType(actualShapeType);
            edges << keyEdge;
            lastDrawnCells_ << keyEdge;
        }
        if (isClosedShape)
        {
            auto keyEdge = vac_->newKeyEdge(interactiveTime(), vertices[verticesCount - 1], vertices[0], 0, w);
            keyEdge->setShapeType(actualShapeType);
            edges << keyEdge;
            lastDrawnCells_ << keyEdge;
        }

        //Draw face
        if (global()->isDrawShapeFaceEnabled() && isClosedShape)
        {
            Cycle cycle(edges);
            auto faceCell = vac_->newKeyFace(cycle);
            faceCell->setColor(global()->faceColor());
            faceCell->setShapeType(actualShapeType);
            lastDrawnCells_ << faceCell->toFaceCell();
        }

        isSuccess = true;
    };

    switch (shapeType) {
    case ShapeType::LINE:
    {
        verticesPoints[0] = QPointF(shapeStartX_, shapeStartY_);
        verticesPoints[1] = QPointF(xScene, yScene);
        processDrawShape(false);
        break;
    }
    case ShapeType::CIRCLE:
    {
        clearLastCells();
        const auto rx = shapeWidth / 2;
        const auto ry = shapeHeight / 2;
        const auto xCenter = leftX + rx;
        const auto yCenter = topY + ry;
        vac_->beginSketchEdge(xCenter + rx, yCenter, w, interactiveTime());
        for (auto i = 0; i <= countAngles; i++)
        {
            auto angle = 2 * M_PI * i / countAngles;
            auto pos = QPointF(rx * cos(angle) + xCenter, ry * sin(angle) + yCenter);
            vac_->continueSketchEdge(pos.x(), pos.y(), w);
        }
        vac_->endSketchEdge(false);
        auto keyEdge = vac_->instantEdges().last();
        keyEdge->setShapeType(shapeType);
        lastDrawnCells_ << keyEdge;

        if (global()->isDrawShapeFaceEnabled())
        {
            vac_->addToSelection(keyEdge, false);
            scene()->createFace(false);
            auto faceCell = vac_->faces().last();
            vac_->addToSelection(faceCell);
            faceCell->setShapeType(shapeType);
            lastDrawnCells_ << faceCell;
        }
        break;
    }
    case ShapeType::TRIANGLE:
    {
        verticesPoints[0] = QPointF(leftX + shapeWidth / 2, topY);
        verticesPoints[1] = QPointF(rightX, bottomY);
        verticesPoints[2] = QPointF(leftX, bottomY);

        processDrawShape();
        break;
    }
    case ShapeType::RECTANGLE:
    {
        verticesPoints[0] = QPointF(leftX, topY);
        verticesPoints[1] = QPointF(rightX, topY);
        verticesPoints[2] = QPointF(rightX, bottomY);
        verticesPoints[3] = QPointF(leftX, bottomY);

        processDrawShape();
        break;
    }
    case ShapeType::POLYGON:
    {
        const auto rx = shapeWidth / 2;
        const auto ry = shapeHeight / 2;
        const auto xCenter = leftX + rx;
        const auto yCenter = topY + ry;
        for (auto i = 0; i < countAngles; i++)
        {
            auto angle = qDegreesToRadians(initialRotation) + 2 * M_PI * i / countAngles;
            auto pos = QPointF(rx * cos(angle) + xCenter, ry * sin(angle) + yCenter);
            verticesPoints[i] = pos;
        }
        processDrawShape();
        break;
    }
    default:
        break;
    }

    lastMousePos_ = currentMousePos;

    if (isSuccess) {
        // Update the selection geometry for display the correct selection size when drawing
        global()->updateSelectedGeometry(leftX, topY, shapeWidth, shapeHeight, 0.0);
    }
}

void View::updateView()
{
    emit allViewsNeedToUpdatePicking();
    updateHoveredObject(mouse_Event_X_, mouse_Event_Y_);
    emit allViewsNeedToUpdate();
}

// Move the selection on delta in millimeters
void View::moveSelected(const double dx_mm, const double dy_mm)
{
    if (vac_->hasSelectedCells()) {
        const auto delta = QPointF(global()->surfaceScaleFactor() * dx_mm, global()->surfaceScaleFactor() * dy_mm);

        vac_->setHoveredCell(vac_->selectedCells().toList().at(0));
        auto bb = vac_->selectedOutlineBoundingBox();
        vac_->prepareDragAndDrop(bb.xMin(), bb.yMin(), interactiveTime());

        // Check if a new position of the shape will be inside the surfaces area
        if (global()->isShapeInSurface(vac_->draggedVertices(), vac_->draggedEdges(), delta)) {
            vac_->performDragAndDrop(bb.xMin() + delta.x(), bb.yMin() + delta.y());
            vac_->completeDragAndDrop();
            vac_->setNoHoveredAllCells();
            updateView();
        }
    }
}

void View::forcedMouseRelease()
{
    QMouseEvent event(QEvent::MouseButtonRelease, lastMousePos_, Qt::LeftButton, Qt::NoButton, Qt::NoModifier);
    qApp->sendEvent(this, &event);
}

void View::drawScene()
{
    if(!mouse_HideCursor_)
    {
        switch(global()->toolMode())
        {
        case Global::SELECT:
            setCursor(Qt::ArrowCursor);
            break;
        case Global::SKETCH:
        case Global::PAINT:
        case Global::SCULPT:
        case Global::DRAW_LINE:
        case Global::DRAW_RECTANGLE:
        case Global::DRAW_CIRCLE:
        case Global::DRAW_TRIANGLE:
        case Global::DRAW_RHOMBUS:
        case Global::DRAW_PENTAGON:
        case Global::DRAW_HEXAGON:
        case Global::DRAW_HEPTAGON:
        case Global::DRAW_OCTAGON:
            setCursor(isMouseInSurface_ ? Qt::CrossCursor : Qt::ArrowCursor);
            break;
        default:
            setCursor(Qt::CrossCursor);
            break;
        }
    }

    // Clear to transparent
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);

    // Note:
    // It is the responsability of view to decide when to call scene()->drawCanvas
    // draw a canvas, and layer backgrounds, and in which order,
    // since this is dependent on onion skinning settings which only
    // view should be aware of

    // Draw canvas
    // XXX Should be replaced by drawCanvas_(scene_->canvas());
    scene_->drawCanvas(viewSettings_);

    // Draw scene
    drawSceneDelegate_(activeTime());
}

void View::drawSceneDelegate_(Time t)
{
    for (int j = 0; j < scene()->numLayers(); ++j)
    {
        Layer * layer = scene()->layer(j);
        if (!layer->isVisible()) {
            continue;
        }
        Background * background = layer->background();
        VectorAnimationComplex::VAC * vac = layer->vac();

        // Draw background
        drawBackground_(background, t.frame());

        // Loop over all onion skins. Draw in this order:
        //   1. onion skins before
        //   2. onion skins after
        //   3. current frame
        //
        // Note 1: For now, we show onions skins for all layers
        //         In the future, by default, we should show onion skins only
        //         for the active layer, and allow user to show them for all
        //         layer via a user option in the onion skin menu.
        //
        // Note 2: Backgrounds are always ignored for onion skinning

        // Draw onion skins
        viewSettings_.setMainDrawing(false);
        if(viewSettings_.onionSkinningIsEnabled())
        {
            // Draw onion skins before
            Time tOnion = t;
            for(int i=0; i<viewSettings_.numOnionSkinsBefore(); ++i)
            {
                tOnion = tOnion - viewSettings_.onionSkinsTimeOffset();
                glTranslated(-viewSettings_.onionSkinsXOffset(),-viewSettings_.onionSkinsYOffset(),0);
            }
            for(int i=0; i<viewSettings_.numOnionSkinsBefore(); ++i)
            {
                vac->draw(tOnion, viewSettings_);
                tOnion = tOnion + viewSettings_.onionSkinsTimeOffset();
                glTranslated(viewSettings_.onionSkinsXOffset(),viewSettings_.onionSkinsYOffset(),0);
            }

            // Draw onion skins after
            tOnion = t;
            for(int i=0; i<viewSettings_.numOnionSkinsAfter(); ++i)
            {
                glTranslated(viewSettings_.onionSkinsXOffset(),viewSettings_.onionSkinsYOffset(),0);
                tOnion = tOnion + viewSettings_.onionSkinsTimeOffset();
                vac->draw(tOnion, viewSettings_);
            }
            for(int i=0; i<viewSettings_.numOnionSkinsAfter(); ++i)
            {
                glTranslated(-viewSettings_.onionSkinsXOffset(),-viewSettings_.onionSkinsYOffset(),0);
            }
        }

        // Draw current frame
        viewSettings_.setMainDrawing(true);
        vac->draw(t, viewSettings_);
    }
}

void View::toggleOutline()
{
    viewSettings_.toggleOutline();
    viewSettingsWidget_->updateWidgetFromSettings();
    update();
}

void View::toggleOutlineOnly()
{
    viewSettings_.toggleOutlineOnly();
    viewSettingsWidget_->updateWidgetFromSettings();
    update();
}

void View::setDisplayMode(ViewSettings::DisplayMode displayMode)
{
    viewSettings_.setDisplayMode(displayMode);
    viewSettingsWidget_->updateWidgetFromSettings();
    update();
}

void View::setOnionSkinningEnabled(bool enabled)
{
    viewSettings_.setOnionSkinningIsEnabled(enabled);
    viewSettingsWidget_->updateWidgetFromSettings();
    update();
}

void View::fitAllInWindow()
{
    // TODO
}

void View::fitSelectionInWindow()
{
    // TODO
}

double View::zoom() const
{
    return camera2D().zoom();
}

// Note: In the future, when rotation of the viewport is allowed,
//       then it should be replaced by:
//           xSceneMin = min(x1, x2, x3, x4);
//           xSceneMax = max(x1, x2, x3, x4);
//           ySceneMin = min(y1, y2, y3, y4);
//           ySceneMax = max(y1, y2, y3, y4);
//       where the (xi,yi)'s are the four corners of the viewport in
//       scene coordinate, which in general will not be axis-aligned

double  View::xSceneMin() const
{
    return - camera2D().x() / zoom();
}

double  View::ySceneMin() const
{
    return - camera2D().y() / zoom();
}

double  View::xSceneMax() const
{
    double x = xSceneMin();
    double w = width();
    double z = zoom();

    return x+w/z;
}

double  View::ySceneMax() const
{
    double x = ySceneMin();
    double w = height();
    double z = zoom();

    return x+w/z;
}

ViewSettings View::viewSettings() const
{
    return viewSettings_;
}

ViewSettingsWidget * View::viewSettingsWidget() const
{
    return viewSettingsWidget_;
}


/***********************************************************
 *              PICKING
 */

void View::drawPick()
{
    Time t = activeTime();
    {
        if(viewSettings_.onionSkinningIsEnabled() && viewSettings_.areOnionSkinsPickable())
        {
            Time tOnion = t;
            for(int i=0; i<viewSettings_.numOnionSkinsBefore(); ++i)
            {
                tOnion = tOnion - viewSettings_.onionSkinsTimeOffset();
                glTranslated(-viewSettings_.onionSkinsXOffset(),-viewSettings_.onionSkinsYOffset(),0);
            }
            for(int i=0; i<viewSettings_.numOnionSkinsBefore(); ++i)
            {
                scene_->drawPick(tOnion, viewSettings_);
                tOnion = tOnion + viewSettings_.onionSkinsTimeOffset();
                glTranslated(viewSettings_.onionSkinsXOffset(),viewSettings_.onionSkinsYOffset(),0);
            }

            tOnion = t;
            for(int i=0; i<viewSettings_.numOnionSkinsAfter(); ++i)
            {
                glTranslated(viewSettings_.onionSkinsXOffset(),viewSettings_.onionSkinsYOffset(),0);
                tOnion = tOnion + viewSettings_.onionSkinsTimeOffset();
                scene_->drawPick(tOnion, viewSettings_);
            }
            for(int i=0; i<viewSettings_.numOnionSkinsAfter(); ++i)
            {
                glTranslated(-viewSettings_.onionSkinsXOffset(),-viewSettings_.onionSkinsYOffset(),0);
            }
        }

        // Draw current frame
        scene_->drawPick(t, viewSettings_);
    }
}

bool View::updateHoveredObject(int x, int y)
{
    // make sure this does NOT redraw the scene, just change its highlighted status.

    if(!pickingIsEnabled_)
        return false;

    // Don't do anything if no picking image
    if(!pickingImg_)
        return false;

    // Find object under the mouse
    Picking::Object old = hoveredObject_;
    if(x<0 || x>=pickingWidth_ || y<0 || y>=pickingHeight_)
    {
        hoveredObject_ = Picking::Object();
    }
    else
    {
        hoveredObject_ = getCloserObject(x, y);
    }

    // Check if it has changed
    bool hasChanged = !(hoveredObject_ == old);

    // If it has, inform the scene of the new highlighted state
    if(hasChanged)
    {
        if(hoveredObject_.isNull())
        {
            scene_->setNoHoveredObject();
        }
        else
        {
            scene_->setHoveredObject(
                activeTime(),
                hoveredObject_.index(),
                hoveredObject_.id());
        }
    }
    else
    {
        // even if it has not changed, inform the scene when nothing's highlighted
        if(hoveredObject_.isNull())
            scene_->setNoHoveredObject();
    }

    return hasChanged;
}

uchar * View::pickingImg(int x, int y)
{
    int k = 4*( (pickingHeight_ - y - 1)*pickingWidth_ + x);
    return &pickingImg_[k];
}

// This method must be very fast. Assumes x and y in range
Picking::Object View::getCloserObject(int x, int y)
{
    // First look directly whether there's an object right at mouse position
    uchar * p = pickingImg(x,y);
    uchar r=p[0], g=p[1], b=p[2];
    if(r!=255 || g!=255 || b!=255)
    {
        return Picking::objectFromRGB(r,g,b);
    }
    else
    {
        // If not, look around in a radius of D pixels
        int D = 3;

        // Clipping
        if(x<D)
            D = x;
        if(y<D)
            D = y;
        int rightBorderDist = pickingWidth_-1-x;
        if(rightBorderDist<D)
            D = rightBorderDist;
        int bottomBorderDist = pickingHeight_-1-y;
        if(bottomBorderDist<D)
            D = bottomBorderDist;


        for(int d=1; d<=D; d++)
        {
            // top row
            for(int varX=x-d; varX<=x+d; varX++)
            {
                p = pickingImg(varX,y-d);
                r=p[0], g=p[1], b=p[2];
                if(r!=255 || g!=255 || b!=255)
                    return Picking::objectFromRGB(r,g,b);
            }
            // bottom row
            for(int varX=x-d; varX<=x+d; varX++)
            {
                p = pickingImg(varX,y+d);
                r=p[0], g=p[1], b=p[2];
                if(r!=255 || g!=255 || b!=255)
                    return Picking::objectFromRGB(r,g,b);
            }
            // left column
            for(int varY=y-d; varY<=y+d; varY++)
            {
                p = pickingImg(x-d,varY);
                r=p[0], g=p[1], b=p[2];
                if(r!=255 || g!=255 || b!=255)
                    return Picking::objectFromRGB(r,g,b);
            }
            // right column
            for(int varY=y-d; varY<=y+d; varY++)
            {
                p = pickingImg(x+d,varY);
                r=p[0], g=p[1], b=p[2];
                if(r!=255 || g!=255 || b!=255)
                    return Picking::objectFromRGB(r,g,b);
            }
        }

        // If still no object found, return a null object
        return Picking::Object();
    }
}

void View::deletePicking()
{
    if(pickingImg_)
    {
        gl_fbo_->glDeleteFramebuffers(1, &fboId_);
        gl_fbo_->glDeleteRenderbuffers(1, &rboId_);
        glDeleteTextures(1, &textureId_);
        hoveredObject_ = Picking::Object();
        delete[] pickingImg_;
        pickingImg_ = 0;
        pickingWidth_ = 0;
        pickingHeight_ = 0;
    }
}

void View::newPicking()
{
    pickingWidth_ = width();
    pickingHeight_ = height();
    pickingImg_ = new uchar[4 * pickingWidth_ * pickingHeight_];

    //  code adapted from http://www.songho.ca/opengl/gl_fbo.html

    // create a texture object
    glGenTextures(1, &textureId_);
    glBindTexture(GL_TEXTURE_2D, textureId_);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE); // automatic mipmap
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, pickingWidth_, pickingHeight_, 0,
                 GL_RGBA, GL_UNSIGNED_BYTE, 0);
    glBindTexture(GL_TEXTURE_2D, 0);

    // create a renderbuffer object to store depth info
    gl_fbo_->glGenRenderbuffers(1, &rboId_);
    gl_fbo_->glBindRenderbuffer(GL_RENDERBUFFER, rboId_);
    gl_fbo_->glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT,
                                   pickingWidth_, pickingHeight_);
    gl_fbo_->glBindRenderbuffer(GL_RENDERBUFFER, 0);

    // create a framebuffer object
    gl_fbo_->glGenFramebuffers(1, &fboId_);
    gl_fbo_->glBindFramebuffer(GL_FRAMEBUFFER, fboId_);

    // attach the texture to FBO color attachment point
    gl_fbo_->glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                                    GL_TEXTURE_2D, textureId_, 0);

    // attach the renderbuffer to depth attachment point
    gl_fbo_->glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                                       GL_RENDERBUFFER, rboId_);

    // check FBO status
    GLenum status = gl_fbo_->glCheckFramebufferStatus(GL_FRAMEBUFFER);
    if(status != GL_FRAMEBUFFER_COMPLETE)
    {
        qDebug() << "ERROR void View::newPicking()"
               << "FBO status != GL_FRAMEBUFFER_COMPLETE";
        return;
    }

    // switch back to window-system-provided framebuffer
    gl_fbo_->glBindFramebuffer(GL_FRAMEBUFFER, defaultFramebufferObject());
}

#include <QElapsedTimer>
#include <QtDebug>

void View::enablePicking()
{
    pickingIsEnabled_ = true;
}

void View::disablePicking()
{
    pickingIsEnabled_ = false;
}

namespace
{

void imageCleanupHandler(void * info)
{
    uchar * img = reinterpret_cast<uchar*>(info);
    delete[] img;
}

}

QImage View::drawToImage(double x, double y, double w, double h, int imgW, int imgH, bool useViewSettings)
{
    return drawToImage(activeTime(), x, y, w, h, imgW, imgH, useViewSettings);
}

QImage View::drawToImage(Time t, double x, double y, double w, double h, int IMG_SIZE_X, int IMG_SIZE_Y, bool useViewSettings)
{
    // Make this widget's rendering context the current OpenGL context
    makeCurrent();


    // ------------ Create multisample FBO --------------------

    GLuint ms_fboId;
    GLuint ms_ColorBufferId;
    GLuint ms_DepthBufferId;
    GLint  ms_samples;

    // Maximum supported samples
    glGetIntegerv(GL_MAX_SAMPLES, &ms_samples);
    // Create FBO
    gl_fbo_->glGenFramebuffers(1, &ms_fboId);
    gl_fbo_->glBindFramebuffer(GL_FRAMEBUFFER, ms_fboId);
    // Create multisample color buffer
    gl_fbo_->glGenRenderbuffers(1, &ms_ColorBufferId);
    gl_fbo_->glBindRenderbuffer(GL_RENDERBUFFER, ms_ColorBufferId);
    gl_fbo_->glRenderbufferStorageMultisample(GL_RENDERBUFFER, ms_samples, GL_RGBA8, IMG_SIZE_X, IMG_SIZE_Y);
    // Create multisample depth buffer
    gl_fbo_->glGenRenderbuffers(1, &ms_DepthBufferId);
    gl_fbo_->glBindRenderbuffer(GL_RENDERBUFFER, ms_DepthBufferId);
    gl_fbo_->glRenderbufferStorageMultisample(GL_RENDERBUFFER, ms_samples, GL_DEPTH_COMPONENT24, IMG_SIZE_X, IMG_SIZE_Y);
    // Attach render buffers to FBO
    gl_fbo_->glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, ms_ColorBufferId);
    gl_fbo_->glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, ms_DepthBufferId);
    // Check FBO status
    GLenum ms_status = gl_fbo_->glCheckFramebufferStatus(GL_FRAMEBUFFER);
    if(ms_status != GL_FRAMEBUFFER_COMPLETE) {
        qDebug() << "Error: FBO ms_status != GL_FRAMEBUFFER_COMPLETE";
        return QImage();
    }


    // ------------ Create standard FBO --------------------

    GLuint fboId;
    GLuint textureId;
    GLuint rboId;

    // Create FBO
    gl_fbo_->glGenFramebuffers(1, &fboId);
    gl_fbo_->glBindFramebuffer(GL_FRAMEBUFFER, fboId);
    // Create color texture
    glGenTextures(1, &textureId);
    glBindTexture(GL_TEXTURE_2D, textureId);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, IMG_SIZE_X, IMG_SIZE_Y, 0,
                 GL_RGBA, GL_UNSIGNED_BYTE, 0);
    glBindTexture(GL_TEXTURE_2D, 0);
    // Create depth buffer
    gl_fbo_->glGenRenderbuffers(1, &rboId);
    gl_fbo_->glBindRenderbuffer(GL_RENDERBUFFER, rboId);
    gl_fbo_->glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, IMG_SIZE_X, IMG_SIZE_Y);
    gl_fbo_->glBindRenderbuffer(GL_RENDERBUFFER, 0);
    // Attach render buffers / textures to FBO
    gl_fbo_->glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textureId, 0);
    gl_fbo_->glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rboId);
    // Check FBO status
    GLenum status = gl_fbo_->glCheckFramebufferStatus(GL_FRAMEBUFFER);
    if(status != GL_FRAMEBUFFER_COMPLETE) {
        qDebug() << "Error: FBO status != GL_FRAMEBUFFER_COMPLETE";
        return QImage();
    }


    // ------------ Render scene to multisample FBO --------------------

    // Bind FBO
    gl_fbo_->glBindFramebuffer(GL_FRAMEBUFFER, ms_fboId);

    // Set viewport
    GLint oldViewport[4];
    glGetIntegerv(GL_VIEWPORT, oldViewport);
    glViewport(0, 0, IMG_SIZE_X, IMG_SIZE_Y);

    // Clear FBO to fully transparent
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Set projection matrix
    // Note: (0,h) and not (h,0) since y-axis is down in VPaint, up in QImage
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, w, 0, h, 0, 1);

    // Set view matrix
    glMatrixMode (GL_MODELVIEW);
    GLWidget_Camera2D camera2d;
    camera2d.setX(-x);
    camera2d.setY(-y);
    camera2d.setZoom(1);
    glLoadMatrixd(camera2d.viewMatrixData());

    // Draw scene
    if (useViewSettings)
    {
        drawSceneDelegate_(t);
    }
    else
    {
        ViewSettings::DisplayMode oldDM = viewSettings_.displayMode();
        viewSettings_.setDisplayMode(ViewSettings::ILLUSTRATION);
        viewSettings_.setMainDrawing(false);
        viewSettings_.setDrawCursor(false);

        for (int j = 0; j < scene()->numLayers(); ++j)
        {
            Layer * layer = scene()->layer(j);
            if (layer->isVisible()) {
                drawBackground_(layer->background(), t.frame());
                layer->vac()->draw(t, viewSettings_);
            }
        }

        viewSettings_.setDrawCursor(true);
        viewSettings_.setDisplayMode(oldDM);
    }

    // Restore viewport
    glViewport(oldViewport[0], oldViewport[1], oldViewport[2], oldViewport[3]);

    // Unbind FBO
    gl_fbo_->glBindFramebuffer(GL_FRAMEBUFFER, defaultFramebufferObject());


    // ------ Blit multisample FBO to standard FBO ---------

    // Bind multisample FBO for reading
    gl_fbo_->glBindFramebuffer(GL_READ_FRAMEBUFFER, ms_fboId);
    // Bind standard FBO for drawing
    gl_fbo_->glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fboId);
    // Blit
    gl_fbo_->glBlitFramebuffer(0, 0, IMG_SIZE_X, IMG_SIZE_Y, 0, 0, IMG_SIZE_X, IMG_SIZE_Y, GL_COLOR_BUFFER_BIT, GL_NEAREST);
    // Unbind FBO
    gl_fbo_->glBindFramebuffer(GL_FRAMEBUFFER, defaultFramebufferObject());


    // ------ Read standard FBO to RAM data ---------

    // Bind standard FBO for reading
    glBindTexture(GL_TEXTURE_2D, textureId);
    // Read
    uchar * img = new uchar[4 * IMG_SIZE_X * IMG_SIZE_Y];
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE,  img);
    // Unbind FBO
    glBindTexture(GL_TEXTURE_2D, 0);


    // ------ Release allocated GPU memory  ---------

    gl_fbo_->glDeleteFramebuffers(1, &ms_fboId);
    gl_fbo_->glDeleteRenderbuffers(1, &ms_ColorBufferId);
    gl_fbo_->glDeleteRenderbuffers(1, &ms_DepthBufferId);
    gl_fbo_->glDeleteFramebuffers(1, &fboId);
    gl_fbo_->glDeleteRenderbuffers(1, &rboId);
    glDeleteTextures(1, &textureId);


    // ------ un-premultiply alpha ---------

    // Once can notice that glBlendFuncSeparate(alpha, 1-alpha, 1, 1-alpha)
    // performs the correct blending function with input:
    //    Frame buffer color as pre-multiplied alpha
    //    Input fragment color as post-multiplied alpha
    // and output:
    //    New frame buffer color as pre-multiplied alpha
    //
    // So by starting with glClearColor(0.0, 0.0, 0.0, 0.0), which is the
    // correct pre-multiplied representation for fully transparent, then
    // by specifying glColor() in post-multiplied alpha, we get the correct
    // blending behaviour and simply have to un-premultiply the value obtained
    // in the frame buffer at the very end

    for(int k=0; k<IMG_SIZE_X*IMG_SIZE_Y; ++k)
    {
        uchar * pixel = &(img[4*k]);
        double a = pixel[3];
        if( 0 < a && a < 255 )
        {
            double s = 255.0 / a;
            pixel[0] = (uchar) (std::min(255.0,std::floor(0.5+s*pixel[0])));
            pixel[1] = (uchar) (std::min(255.0,std::floor(0.5+s*pixel[1])));
            pixel[2] = (uchar) (std::min(255.0,std::floor(0.5+s*pixel[2])));
        }
    }


    // ------ Convert to Qimage ---------

    // Create cleanup info to delete[] img when appropriate
    QImageCleanupFunction cleanupFunction = &imageCleanupHandler;
    void * cleanupInfo = reinterpret_cast<void*>(img);

    // Create QImage
    QImage res(img, IMG_SIZE_X, IMG_SIZE_Y, QImage::Format_RGBA8888, cleanupFunction, cleanupInfo);

    // Return QImage
    return res;
}

void View::updatePicking()
{
    // Remove previously highlighted object
    hoveredObject_ = Picking::Object();

    if(!pickingIsEnabled_)
        return;

    // Make this widget's rendering context the current OpenGL context
    makeCurrent();

    // get the viewport size, allocate memory if necessary
    if( !(width()>0) || !(height()>0))
    {
        deletePicking();
        return;
    }
    else if(
        pickingImg_
        && (pickingWidth_ == width())
        && (pickingHeight_ == height()))
    {
        // necessary objects already created: do nothing
    }
    else
    {
        deletePicking();
        newPicking();
    }

    // set rendering destination to FBO
    gl_fbo_->glBindFramebuffer(GL_FRAMEBUFFER, fboId_);

    // clear buffers
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Should we setup other things? (e.g., disabling antialiasing)
    // Seems to work as is. If issues, check GLWidget::initilizeGL()

    // Set viewport
    GLint oldViewport[4];
    glGetIntegerv(GL_VIEWPORT, oldViewport);
    glViewport(0, 0, pickingWidth_, pickingHeight_);

    // Setup camera position and orientation
    setCameraPositionAndOrientation();

    // draw the picking
    drawPick();

    // Restore viewport
    glViewport(oldViewport[0], oldViewport[1], oldViewport[2], oldViewport[3]);

    // unbind FBO
    gl_fbo_->glBindFramebuffer(GL_FRAMEBUFFER, defaultFramebufferObject());

    // extract the texture info from GPU to RAM: EXPENSIVE + MAY CAUSE OPENGL STALL
    glBindTexture(GL_TEXTURE_2D, textureId_);
    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, pickingImg_);
    glBindTexture(GL_TEXTURE_2D, 0);

    // Update highlighted object
    if(underMouse())
    {
        updateHoveredObject(mouse_Event_X_, mouse_Event_Y_);
    }
}
