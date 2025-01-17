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

#ifndef TIMELINE_H
#define TIMELINE_H

#include <QWidget>
#include <QDialog>
#include <QComboBox>
#include <QElapsedTimer>
#include <QList>
#include <QSet>
#include <QColor>
#include "TimeDef.h"
#include "VAC/vpaint_global.h"

class QPushButton;
class QSpinBox;
class QTimer;
class QHBoxLayout;
class QCheckBox;
class View;
namespace VPaint
{
class Scene;
}
class XmlStreamWriter;
class XmlStreamReader;
class QAction;

// Paint the top timeline bar.
// It's a friend of Timeline: can access all members of it
class Timeline;
class Q_VPAINT_EXPORT Timeline_HBar: public QWidget
{
public:
    Timeline_HBar(Timeline * w);

protected:
    virtual void paintEvent (QPaintEvent * event);
    virtual void mouseMoveEvent (QMouseEvent * event);
    virtual void mousePressEvent (QMouseEvent * event);
    virtual void mouseReleaseEvent (QMouseEvent * event);
    virtual void leaveEvent (QEvent * event);
    
private:
    Timeline * w_;
    bool isScrolling_;
    int scrollingInitialX_;
    int scrollingInitialFrame_;
    int scrollingInitialOffset_;

    bool hasHighlightedFrame_;
    int highlightedFrame_;

    QList<QColor> colors_;
};

class PlaybackSettings
{
public:
    PlaybackSettings();
    void setDefaultValues(); // Restore default values / Reset

    enum PlayMode {
        NORMAL = 0,
        LOOP,
        BOUNCE
    };
    static QString playModeToString(PlayMode mode);
    PlayMode stringToPlayMode(const QString & str);

    int firstFrame() const;
    int lastFrame() const;
    int fps() const;
    PlayMode playMode() const;
    bool subframeInbetweening() const;

    void setFirstFrame(int f);
    void setLastFrame(int f);
    void setFps(int n) ;
    void setPlayMode(PlayMode mode);
    void setSubframeInbetweening(bool b);

    void read(XmlStreamReader & xml);
    void write(XmlStreamWriter & xml) const;

private:
    int firstFrame_;
    int lastFrame_;
    int fps_;
    PlayMode playMode_;
    bool subframeInbetweening_;
};

class Q_VPAINT_EXPORT PlaybackSettingsDialog: public QDialog
{
public:
    PlaybackSettingsDialog(const PlaybackSettings & settings = PlaybackSettings());

    PlaybackSettings playbackSettings() const;
    void setPlaybackSettings(const PlaybackSettings & settings);

private:
    mutable PlaybackSettings settings_;

    QSpinBox * fpsSpinBox_;
    QCheckBox * subframeCheckBox_;
    QComboBox * playModeSpinBox_;
};


class Q_VPAINT_EXPORT Timeline : public QWidget
{
    Q_OBJECT

public:
    Timeline(VPaint::Scene * scene, QWidget *parent = 0);
    ~Timeline();

    void read(XmlStreamReader & xml);
    void write(XmlStreamWriter & xml) const;

    // Set info about the selected cell
    void setSelectionType(int type); // 0 for no selection; 1 for at least one key cells; 2 for only inbetween cells
    void setT(double t);
    void setT1(double t1);
    void setT2(double t2);

    // Set which View's time this timeline control
    void addView(View * view);
    void removeView(View * view);

    // Get playback settings
    int firstFrame() const;
    int lastFrame() const;
    int fps() const;
    PlaybackSettings::PlayMode playMode() const;
    bool subframeInbetweening() const;

    // Current state
    bool isPlaying() const;
    QSet<View*> playedViews() const;

    // Visualization
    int firstVisibleFrame() const;
    int lastVisibleFrame() const;

    QAction * actionGoToFirstFrame() const;
    QAction * actionGoToPreviousFrame() const;
    QAction * actionPlayPause() const;
    QAction * actionGoToNextFrame() const;
    QAction * actionGoToLastFrame() const;

public slots:
    void play();
    void pause();
    void playPause();

    void openPlaybackSettingsDialog();

    void goToFirstFrame();
    void goToPreviousFrame();
    void goToNextFrame();
    void goToLastFrame();

    void setFirstFrame(int firstFrame);
    void setLastFrame(int lastFrame);
    void setFps(int fps);
    void realTimePlayingChanged();

private slots:
    void goToFirstFrame(View * view);
    void goToPreviousFrame(View * view);
    void goToNextFrame(View * view);
    void goToLastFrame(View * view);
    void goToFrame(View * view, int frame);
    void goToFrame(View * view, double frame);

    void timerTimeout();
    void roundPlayedViews();

signals:
    void timeChanged();
    void playingWindowChanged();

protected:
    void paintEvent(QPaintEvent * event);

private:
    // Selected cell info
    int selectionType_; // 0 for no selection; 1 for at least one key cells; 2 for only inbetween cells
    double t_;
    double t1_;
    double t2_;

    // Linked scene
    VPaint::Scene * scene_;

    // The views whose times are controlled by this timeline
    QList<View*> views_;
    QSet<View*> playedViews_;
    
    // Delegate timeline painting and mouse events handling
    friend class Timeline_HBar;
    Timeline_HBar * hbar_;

    // Time control
    QTimer * timer_;
    QElapsedTimer elapsedTimer_;

    // Actions
    QAction * actionGoToFirstFrame_;
    QAction * actionGoToPreviousFrame_;
    QAction * actionPlayPause_;
    QAction * actionGoToNextFrame_;
    QAction * actionGoToLastFrame_;

    // playing properties
    QPushButton * firstFrameButton_;
    QPushButton * previousKeyFrameButton_;
    QPushButton * previousFrameButton_;
    QPushButton * playPauseButton_;
    QPushButton * nextFrameButton_;
    QPushButton * nextKeyFrameButton_;
    QPushButton * lastFrameButton_;
    bool isPlaying_;
    
    // Settings
    PlaybackSettings settings_;

    // Playing window
    QSpinBox * firstFrameSpinBox_;
    QSpinBox * lastFrameSpinBox_;

    // Implementation detail of bounce playback mode
    bool playingDirection_;

    // visible window
    int firstVisibleFrame_;
    int lastVisibleFrame_;
    int totalPixelOffset_;
};

    
    
#endif
