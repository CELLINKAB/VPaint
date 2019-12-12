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

#include "SvgImportDialog.h"

#include <QDialogButtonBox>
#include <QLabel>
#include <QVBoxLayout>

SvgImportDialog::SvgImportDialog(QWidget * parent) :
    QDialog(parent)
{
    setWindowTitle(tr("SVG Import"));
    setMinimumSize(300, 200);

    // Warning label.
    //
    // Note: Qt RichText lacks good support for specifying margin/padding of
    // lists (ul, ol). See: https://stackoverflow.com/questions/24000206/
    //
    // Also, Qt doesn't seem to implement style inheritance, that is, child
    // elements don't inherit the style of parent element, at least for
    // margin/padding (seems to work for font attributes, though)
    //
    QLabel * warning = new QLabel(tr(
        "<p style=\"margin:0;padding:0\">Warning! This importer is BETA and does not support:</p>"
        "<ul style=\"-qt-list-indent:0;margin:5px 0px 10px 15px;padding:0;\">"
          "<li>Gradients, patterns, markers, and dashes</li>"
          "<li>Stroke caps and joins other than 'round'</li>"
          "<li>Masking and clipping</li>"
          "<li>Blur and other filters</li>"
          "<li>Clone instancing using symbols</li>"
          "<li>Embedded or linked images</li>"
          "<li>Units other than 'px'</li>"
          "<li>CSS stylesheets</li>"
        "</lu>"));
    warning->setWordWrap(true);
    warning->setTextFormat(Qt::RichText);

    // Dialog button box
    QDialogButtonBox * buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok);
    connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
    connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));

    // Layout
    QVBoxLayout * layout = new QVBoxLayout();
    layout->addWidget(warning);
    layout->addStretch();
    layout->addWidget(buttonBox);
    setLayout(layout);
}
