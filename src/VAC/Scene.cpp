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

#include "Timeline.h"
#include "Scene.h"
#include "SceneObject.h"
#include <QKeyEvent>

#include "VectorAnimationComplex/VAC.h"
#include "VectorAnimationComplex/InbetweenFace.h"
#include "Background/Background.h"

#include <QtDebug>

#include "XmlStreamReader.h"
#include "XmlStreamWriter.h"

#include "OpenGL.h"
#include "Global.h"

#include "SaveAndLoad.h"

#include <QToolBar>

#include "Layer.h"
namespace VPaint {
Scene::Scene() :
    activeLayerIndex_(-1),
    left_(0),
    top_(0),
    width_(1280),
    height_(720),
    isUseConsistentLayerHeight_(true),
    consistentLayerHeight_(-1.0),
    firstLayerHeightPercents_(-1.0),
    loadedSurfaceType_(SurfaceType::None)
{
    indexHovered_ = -1;
}

Scene * Scene::createDefaultScene()
{
    Scene * res = new Scene();
    auto layer = res->createLayer(tr("Layer 1"));
    layer->background()->setOpacity(1.0);
    return res;
}

double Scene::left() const
{
    return left_;
}

double Scene::top() const
{
    return top_;
}

double Scene::width() const
{
    return width_;
}

double Scene::height() const
{
    return height_;
}

void Scene::setLeft(double x)
{
    left_ = x;
    emitChanged();
}

void Scene::setTop(double y)
{
    top_ = y;
    emitChanged();
}

void Scene::setWidth(double w)
{
    width_ = w;
    emitChanged();
}

void Scene::setHeight(double h)
{
    height_ = h;
    emitChanged();
}

void Scene::setCanvasDefaultValues()
{
    left_ = 0;
    top_ = 0;
    width_ = 1280;
    height_ = 720;

    // Don't emit changed on purpose
}

void Scene::setInfillDensityForSelectedCells(int density)
{
    if(const auto layer = activeLayer())
    {
        layer->vac()->setInfillDensityForSelectedCells(density);
    }
}

void Scene::setInfillPatternForSelectedCells(InfillPattern::Pattern pattern)
{
    if(const auto layer = activeLayer())
    {
        layer->vac()->setInfillPatternForSelectedCells(pattern);
    }
}

void Scene::copyFrom(Scene * other)
{
    // XXX We should also copy canvas properties

    // Block signals
    blockSignals(true);

    // Reset to default
    clear(true);

    // Copy layers
    for(Layer * layer: qAsConst(other->layers_))
        addLayer_(layer->clone(), true);
    activeLayerIndex_ = other->activeLayerIndex_;

    // Reset hovered
    indexHovered_ = -1;

    // Unblock signals
    blockSignals(false);

    // Emit signals
    emit needUpdatePicking();
    emitChanged();
    emit selectionChanged();
    emit layerAttributesChanged();
}

void Scene::clear(bool silent)
{
    // Block signals
    blockSignals(true);

    // Delete all layers
    for(Layer * layer: qAsConst(layers_))
        delete layer;
    layers_.clear();
    activeLayerIndex_ = -1;

    // XXX Shouldn't this clear left/top/width/height too?

    // Unblock signals
    blockSignals(false);

    // Emit signals
    if(!silent)
    {
        emitChanged();
        emit needUpdatePicking();
        emit selectionChanged();
        emit layerAttributesChanged();
    }
}

Scene::~Scene()
{
    clear(true);
}

// ----------------------- Save and Load -------------------------

void Scene::save(QTextStream & out)
{
    Q_UNUSED(out);

    // XXX Deprecated
    /*
    out << Save::newField("SceneObjects");
    out << "\n" << Save::indent() << "[";
    Save::incrIndent();
    for(Layer * layer: layers_)
    {
        out << Save::openCurlyBrackets();
        layer->vac()->save(out);
        out << Save::closeCurlyBrackets();
    }
    Save::decrIndent();
    out << "\n" << Save::indent() << "]";
    */
}

void Scene::exportSVG(Time t, QTextStream & out)
{
    // Export Layers
    for(Layer * layer: qAsConst(layers_))
    {
        layer->background()->exportSVG(
            t.frame(), out, left(), top(), width(), height());
        layer->exportSVG(t, out);
    }
}

void Scene::read(QTextStream & in)
{
    Q_UNUSED(in);

    // XXX Deprecated
    /*
    clear(true);
    
    QString field;
    QString bracket;
    field = Read::field(in);
    Read::skipBracket(in); // [
    while(Read::string(in) == "{")
    {
        addLayer(SceneObject::read(in), true);
        Read::skipBracket(in); // }
    }
    // if here, last read string == ]


    VectorAnimationComplex::VAC * vac = activeLayer();
    if(vac)
    {
        connect(vac,SIGNAL(selectionChanged()),this,SIGNAL(selectionChanged()));
    }

    emit changed();
    emit needUpdatePicking();
    emit selectionChanged();
    */
}

void Scene::writeAllLayers(XmlStreamWriter & xml) const
{
    for(Layer * layer: qAsConst(layers_))
    {
        xml.writeStartElement("layer");
        layer->write(xml);
        xml.writeEndElement();
    }
}

void Scene::readOneLayer(XmlStreamReader & xml)
{
    // Precondition: XML element "layer" just opened

    // XXX Remove this
    blockSignals(true);

    Layer * layer = new Layer();
    layer->read(xml);
    addLayer_(layer, true);

    // XXX Remove these: only emit once after finishing to read the file
    blockSignals(false);
    emit needUpdatePicking();
    emit changed();
    emit selectionChanged();
    emit layerAttributesChanged();
}

void Scene::readCanvas(XmlStreamReader & xml)
{
    setCanvasDefaultValues();

    // Canvas
    if(xml.attributes().hasAttribute("position"))
    {
        QString stringPos = xml.attributes().value("position").toString();
        QStringList list = stringPos.split(" ");
        setLeft(list[0].toDouble());
        setTop(list[1].toDouble());
    }

    if(xml.attributes().hasAttribute("size"))
    {
        QString stringsize = xml.attributes().value("size").toString();
        QStringList list = stringsize.split(" ");
        setWidth(list[0].toDouble());
        setHeight(list[1].toDouble());
    }

    if (xml.attributes().hasAttribute("useConsistentLayerHeight"))
    {
        isUseConsistentLayerHeight_ = xml.attributes().value("useConsistentLayerHeight").toInt();
    }

    if (xml.attributes().hasAttribute("consistentLayerHeight"))
    {
        consistentLayerHeight_ = xml.attributes().value("consistentLayerHeight").toDouble();
    }

    if (xml.attributes().hasAttribute("firstLayerHeightPercents"))
    {
        firstLayerHeightPercents_ = xml.attributes().value("firstLayerHeightPercents").toDouble();
    }

    if (xml.attributes().hasAttribute("surfaceType"))
    {
        loadedSurfaceType_ = static_cast<SurfaceType>(xml.attributes().value("surfaceType").toInt());
    }

    if (xml.attributes().hasAttribute("sceneWidth"))
    {
        setWidth(xml.attributes().value("sceneWidth").toInt());
    }

    if (xml.attributes().hasAttribute("sceneHeight"))
    {
        setHeight(xml.attributes().value("sceneHeight").toInt());
    }


    xml.skipCurrentElement();
}

void Scene::writeCanvas(XmlStreamWriter & xml) const
{
    xml.writeAttribute("position", QString().setNum(left()) + " " + QString().setNum(top()));
    xml.writeAttribute("size", QString().setNum(width()) + " " + QString().setNum(height()));
    xml.writeAttribute("useConsistentLayerHeight", QString().setNum(isUseConsistentLayerHeight()));
    xml.writeAttribute("consistentLayerHeight", QString().setNum(consistentLayerHeight()));
    xml.writeAttribute("firstLayerHeightPercents", QString().setNum(firstLayerHeightPercents()));
    xml.writeAttribute("surfaceType", QString().setNum(static_cast<int>(global()->surfaceType())));
    xml.writeAttribute("sceneWidth", QString().setNum(width()));
    xml.writeAttribute("sceneHeight", QString().setNum(height()));
}

void Scene::relativeRemap(const QDir & oldDir, const QDir & newDir)
{
    for(Layer * layer: qAsConst(layers_))
    {
        layer->background()->relativeRemap(oldDir, newDir);
    }
}

// ----------------------- Drawing the scene -------------------------

// XXX Refactor this: move it to View. Even better, have a Canvas and
// CanvasRenderer class
void Scene::drawCanvas(ViewSettings & /*viewSettings*/)
{
    double x = left();
    double y = top();
    double w = width();
    double h = height();

    if(global()->showCanvas())
    {
        // Out-of-canvas background color
#ifdef CELLINK_VPAINT_STYLE
        glClearColor(0.941f, 0.945f, 0.949f, 1.0f);
#else
        glClearColor(0.8f, 0.8f, 0.8f, 1.0f);
#endif
        glClear(GL_COLOR_BUFFER_BIT);

        // Canvas border
//        glColor4d(0.0, 0.0, 0.0, 1.0);
//        glLineWidth(3.0);
//        glBegin(GL_LINE_LOOP);
//        {
//            glVertex2d(x,y);
//            glVertex2d(x+w,y);
//            glVertex2d(x+w,y+h);
//            glVertex2d(x,y+h);
//        }
//        glEnd();

        // Canvas color
        glColor4d(1.0, 1.0, 1.0, 1.0);
        glBegin(GL_QUADS);
        {
            glVertex2d(x,y);
            glVertex2d(x+w,y);
            glVertex2d(x+w,y+h);
            glVertex2d(x,y+h);
        }
        glEnd();
    }
    else
    {
        // Canvas color
        glColor4d(1.0, 1.0, 1.0, 1.0);
        glClear(GL_COLOR_BUFFER_BIT);
    }
}

void Scene::draw(Time time, ViewSettings & viewSettings)
{
    // Draw layers
    for(Layer * layer: qAsConst(layers_))
    {
        layer->draw(time, viewSettings);
    }
}

void Scene::drawPick(Time time, ViewSettings & viewSettings)
{
    // Find which layer to pick
    Layer * layer = activeLayer();
    int index = -1;
    for(int i=0; i<layers_.size(); i++)
    {
        if (layers_[i] == layer)
        {
            index = i;
            break;
        }
    }

    // Draw picking information
    if (index >= 0)
    {
        Picking::setIndex(index);
        layer->drawPick(time, viewSettings);
    }
}


// ---------------- Highlighting and Selecting -----------------------
    
// No need to emit changed() or needUpdatePicking() here, since highlighting and selecting
// is trigerred by View or View3D, and hence they can decide themselves what do they
// need to update

void Scene::setHoveredObject(Time time, int index, int id)
{
    setNoHoveredObject();
    indexHovered_ = index >= 0 && index < numLayers() ? index : -1;
    if(indexHovered_ != -1)
    {
        layers_[index]->setHoveredObject(time, id);
        layers_[index]->vac()->hoveverShape();
    }
}

void Scene::setNoHoveredObject()
{
    indexHovered_ = indexHovered_ >= 0 && indexHovered_ < numLayers() ? indexHovered_ : -1;
    if(indexHovered_ != -1)
    {
        layers_[indexHovered_]->setNoHoveredObject();
        indexHovered_ = -1;
    }
}

void Scene::select(Time time, int index, int id)
{
    if (index >= 0 && index < numLayers()) {
        layers_[index]->select(time, id);
    }
}

void Scene::deselect(Time time, int index, int id)
{
    if (index >= 0 && index < numLayers()) {
        layers_[index]->deselect(time, id);
    }
}

void Scene::toggle(Time time, int index, int id)
{
    if (index >= 0 && index < numLayers()) {
        layers_[index]->toggle(time, id);
    }
}

void Scene::deselectAll(Time time)
{
    for(Layer * layer: qAsConst(layers_))
        layer->deselectAll(time);
}

void Scene::deselectAll()
{
    for(Layer * layer: qAsConst(layers_))
        layer->deselectAll();
}

void Scene::invertSelection()
{
    for(Layer * layer: qAsConst(layers_))
        layer->invertSelection();
}

void Scene::selectAllInFrame()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->selectAllAtTime(global()->activeTime());
    }
}

void Scene::selectAllInAnimation()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->selectAll();
    }
}

void Scene::selectConnected()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->selectConnected();
    }
}

void Scene::selectClosure()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->selectClosure();
    }
}

void Scene::selectVertices()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->selectVertices();
    }
}

void Scene::selectEdges()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->selectEdges();
    }
}

void Scene::selectFaces()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->selectFaces();
    }
}

void Scene::deselectVertices()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->deselectVertices();
    }
}

void Scene::deselectEdges()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->deselectEdges();
    }
}

void Scene::deselectFaces()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->deselectFaces();
    }
}

void Scene::selectKeyCells()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->selectKeyCells();
    }
}


void Scene::selectInbetweenCells()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->selectInbetweenCells();
    }
}


void Scene::deselectKeyCells()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->deselectKeyCells();
    }
}


void Scene::deselectInbetweenCells()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->deselectInbetweenCells();
    }
}

void Scene::selectKeyVertices()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->selectKeyVertices();
    }
}

void Scene::selectKeyEdges()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->selectKeyEdges();
    }
}

void Scene::selectKeyFaces()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->selectKeyFaces();
    }
}

void Scene::deselectKeyVertices()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->deselectKeyVertices();
    }
}

void Scene::deselectKeyEdges()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->deselectKeyEdges();
    }
}

void Scene::deselectKeyFaces()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->deselectKeyFaces();
    }
}

void Scene::selectInbetweenVertices()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->selectInbetweenVertices();
    }
}

void Scene::selectInbetweenEdges()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->selectInbetweenEdges();
    }
}

void Scene::selectInbetweenFaces()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->selectInbetweenFaces();
    }
}

void Scene::deselectInbetweenVertices()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->deselectInbetweenVertices();
    }
}

void Scene::deselectInbetweenEdges()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->deselectInbetweenEdges();
    }
}

void Scene::deselectInbetweenFaces()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->deselectInbetweenFaces();
    }
}

void Scene::keyPressEvent(QKeyEvent *event)
{
    event->ignore();
}

void Scene::keyReleaseEvent(QKeyEvent *event)
{
    event->ignore();
}

void Scene::addLayer_(Layer * layer, bool silent)
{
    layers_ << layer;
    if (activeLayerIndex_ < 0) {
        activeLayerIndex_ = 0;
    }

    connect(layer, SIGNAL(changed()), this, SIGNAL(changed()));
    connect(layer, SIGNAL(checkpoint()), this, SIGNAL(checkpoint()));
    connect(layer, SIGNAL(needUpdatePicking()), this, SIGNAL(needUpdatePicking()));
    connect(layer, SIGNAL(selectionChanged()), this, SIGNAL(selectionChanged()));
    connect(layer, SIGNAL(layerAttributesChanged()), this, SIGNAL(layerAttributesChanged()));

    if(!silent)
    {
        emitChanged();
        emit needUpdatePicking();
        emit layerAttributesChanged();
    }
}

void Scene::populateToolBar(QToolBar * toolBar)
{
    VectorAnimationComplex::VAC::populateToolBar(toolBar, this);
}

QList<ShapeType>Scene::getActiveLayerShapesType()
{
    QList<ShapeType> shapesType;
    Layer * layer = activeLayer();
    if(layer)
    {
        shapesType =  layer->vac()->getAllShapesType();
    }
    return shapesType;
}

void Scene::setUseConsistentLayerHeight(bool useConsistentHeight)
{
    isUseConsistentLayerHeight_ = useConsistentHeight;
}

void Scene::setConsistentLayerHeight(const qreal height)
{
    consistentLayerHeight_ = height;
}

void Scene::setFirstLayerHeightPercents(const qreal percents)
{
    if (percents >= 0 && percents <= 100) {
        firstLayerHeightPercents_ = percents;
    }
}

SurfaceType Scene::getLoadedSurfaceType()
{
    auto surfaceType = loadedSurfaceType_;
    loadedSurfaceType_ = SurfaceType::None;
    return surfaceType;
}

bool Scene::hasSceneContent() const
{
    for (auto layer : layers_)
    {
        if (layer && layer->vac()->instantVertices().count() > 0) {
            return true;
        }
    }
    return false;
}

bool Scene::isAllShapesInSurface() const
{
    for (auto layer : layers_)
    {
        if (layer && !global()->isShapeInSurface(layer->vac()->instantVertices(), layer->vac()->instantEdges())) {
            return false;
        }
    }
    return true;
}

void Scene::deleteSelectedCells()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        QList<ShapeType> shapesType = layer->vac()->getSelectedShapeType();
        layer->vac()->deleteSelectedCells();
        for(auto shapeType : shapesType)
        if(shapeType != ShapeType::NONE)
        {
            emitShapeDelete();
        }
    }
}

void Scene::test()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->test();
    }
}

void Scene::smartDelete()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        //TOD - multiselection
        QList<ShapeType> shapesType = layer->vac()->getSelectedShapeType();
        layer->vac()->smartDelete();
        for(auto shapeType : shapesType)
        if(shapeType != ShapeType::NONE)
        {
            emitShapeDelete();
        }
    }
}

int Scene::numLayers() const
{
    return layers_.size();
}

Layer* Scene::layer(int i) const
{
    if (0 <= i && i < numLayers())
    {
        return layers_[i];
    }
    else
    {
        return nullptr;
    }
}

void Scene::setActiveLayer(int i, bool needEmitCheckpoint)
{
    if(i != activeLayerIndex_ && 0 <= i && i < numLayers())
    {
        deselectAll();

        activeLayerIndex_ = i;

        emitChanged();
        emit needUpdatePicking();
        emit layerAttributesChanged();
        if (needEmitCheckpoint) {
            emitCheckpoint();
        }
    }
}

Layer * Scene::activeLayer() const
{
    return layer(activeLayerIndex_);
}

int Scene::activeLayerIndex() const
{
    return activeLayerIndex_;
}

VectorAnimationComplex::VAC * Scene::activeVAC() const
{
    Layer * layer = activeLayer();
    return layer ? layer->vac() : nullptr;
}

Background * Scene::activeBackground() const
{
    Layer * layer = activeLayer();
    return layer ? layer->background() : nullptr;
}

Layer * Scene::createLayer()
{
    return createLayer(tr("Layer %1").arg(numLayers() + 1));
}

void Scene::addLayer(Layer * layer , bool setActiveOnTop)
{
    deselectAll();
    addLayer_(layer, true);

    // Move above active layer, or keep last if no active layer
    int newActiveLayerIndex = layers_.size() - 1;
    if (!setActiveOnTop && 0 <= activeLayerIndex_ && activeLayerIndex_ < newActiveLayerIndex - 1)
    {
        newActiveLayerIndex = activeLayerIndex_ + 1;
        for (int i = layers_.size() - 1; i > newActiveLayerIndex; --i)
        {
            layers_[i] = layers_[i-1];
        }
        layers_[newActiveLayerIndex] = layer;
    }

    // Set as active
    activeLayerIndex_ = newActiveLayerIndex;

    // Emit signals
    emitChanged();
    emit needUpdatePicking();
    emit layerAttributesChanged();
    if (layers_.size() > 1) {
        emitCheckpoint();
    }
}

Layer * Scene::createLayer(const QString & name)
{
    // Create new layer, add it on top for now
    Layer * layer = new Layer(name);
    addLayer(layer);
    return layer;
}

void Scene::moveActiveLayerUp()
{
    int i = activeLayerIndex_;
    if(0 <= i && i < numLayers() - 1)
    {
        // Swap out layers
        int j = i + 1;
        std::swap(layers_[i], layers_[j]);

        // Set new active index.
        // Note: we don't call setActiveLayer(j) to avoid calling deselectAll().
        activeLayerIndex_ = j;

        // Emit signals. Note: it may be tempting to think that updating
        // picking is unnecessary (because the pickable cells haven't moved),
        // however the picking image data contains the layer index which has
        // changed, and therefore this picking image needs to be re-rendered.
        emitChanged();
        emit needUpdatePicking();
        emit layerAttributesChanged();
    }
}

void Scene::moveActiveLayerDown()
{
    int i = activeLayerIndex_;
    if(0 < i && i < numLayers())
    {
        // Swap out layers
        int j = i - 1;
        std::swap(layers_[i], layers_[j]);

        // Set new active index.
        // Note: we don't call setActiveLayer(j) to avoid calling deselectAll().
        activeLayerIndex_ = j;

        // Emit signals. Note: it may be tempting to think that updating
        // picking is unnecessary (because the pickable cells haven't moved),
        // however the picking image data contains the layer index which has
        // changed, and therefore this picking image needs to be re-rendered.
        emitChanged();
        emit needUpdatePicking();
        emit layerAttributesChanged();
    }
}

void Scene::destroyActiveLayer()
{
    destroyLayer(activeLayerIndex_);
}

void Scene::destroyLayer(const int index)
{
    if (0 <= index && index < numLayers())
    {
        deselectAll();

        Layer* toBeDestroyedLayer = layers_[index];
        layers_.removeAt(index);

        // Set as active the layer below, unless it was the bottom-most layer
        // or the only layer in the scene.
        if (numLayers() == 0)
        {
            // no more layers
            activeLayerIndex_ = -1;
        }
        else if (activeLayerIndex_ == 0)
        {
            // was the bottom-most layer
            activeLayerIndex_ = 0;
            activeLayer()->background()->setOpacity(1.0);
        }
        else if (index <= activeLayerIndex_)
        {
            // had layer below
            --activeLayerIndex_;
        }

        delete toBeDestroyedLayer;

        emitChanged();
        emit needUpdatePicking();
        emit layerAttributesChanged();
        emitCheckpoint();
    }
}

VectorAnimationComplex::InbetweenFace * Scene::createInbetweenFace()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        return layer->vac()->newInbetweenFace(
                    QList<VectorAnimationComplex::AnimatedCycle>(),
                    QSet<VectorAnimationComplex::KeyFace*>(),
                    QSet<VectorAnimationComplex::KeyFace*>());
    }
    else
    {
        return nullptr;
    }
}

void Scene::cut(VectorAnimationComplex::VAC* & clipboard)
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->cut(clipboard);
    }
}

void Scene::copy(VectorAnimationComplex::VAC* & clipboard)
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->copy(clipboard);
    }
}

void Scene::paste(VectorAnimationComplex::VAC* & clipboard, bool isMousePaste)
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->paste(clipboard, isMousePaste);
    }
}

void Scene::motionPaste(VectorAnimationComplex::VAC* & clipboard)
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->motionPaste(clipboard);
    }
}

void Scene::createFace(bool emitCheckpoint)
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->createFace(emitCheckpoint);
    }
}

void Scene::addCyclesToFace()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->addCyclesToFace();
    }
}

void Scene::removeCyclesFromFace()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->removeCyclesFromFace();
    }
}

void Scene::changeColor()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->changeColor();
    }
}
void Scene::raise()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->raise();
    }
}

void Scene::lower()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->lower();
    }
}

void Scene::raiseToTop()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->raiseToTop();
    }
}

void Scene::lowerToBottom()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->lowerToBottom();
    }
}

void Scene::altRaise()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->altRaise();
    }
}

void Scene::altLower()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->altLower();
    }
}

void Scene::altRaiseToTop()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->altRaiseToTop();
    }
}

void Scene::altLowerToBottom()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->altLowerToBottom();
    }
}

void Scene::changeEdgeWidth()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->changeEdgeWidth();
    }
}

void Scene::glue()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->glue();
    }
}

void Scene::unglue()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->unglue();
    }
}

void Scene::uncut()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->uncut();
    }
}

void Scene::inbetweenSelection()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->inbetweenSelection();
    }
}

void Scene::keyframeSelection()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->keyframeSelection();
    }
}

void Scene::resetCellsToConsiderForCutting()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->resetCellsToConsiderForCutting();
    }
}

void Scene::updateCellsToConsiderForCutting()
{
    Layer * layer = activeLayer();
    if(layer)
    {
        layer->vac()->updateCellsToConsiderForCutting();
    }
}
}
