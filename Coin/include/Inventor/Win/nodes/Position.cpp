// src\Inventor\Win\nodes\Position.cpp.  Generated from Position.cpp.in by configure.

/**************************************************************************\
 *
 *  This file is part of the Coin 3D visualization library.
 *  Copyright (C) by Kongsberg Oil & Gas Technologies.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  ("GPL") version 2 as published by the Free Software Foundation.
 *  See the file LICENSE.GPL at the root directory of this source
 *  distribution for additional information about the GNU GPL.
 *
 *  For using Coin with software that can not be combined with the GNU
 *  GPL, and for taking advantage of the additional benefits of our
 *  support services, please contact Kongsberg Oil & Gas Technologies
 *  about acquiring a Coin Professional Edition License.
 *
 *  See http://www.coin3d.org/ for more information.
 *
 *  Kongsberg Oil & Gas Technologies, Bygdoy Alle 5, 0257 Oslo, NORWAY.
 *  http://www.sim.no/  sales@sim.no  coin-support@coin3d.org
 *
\**************************************************************************/

#include <assert.h>

#include <Inventor/errors/SoDebugError.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoPickAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/actions/SoGetMatrixAction.h>
#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/SoPath.h>

#include <Inventor/Win/nodes/SoGuiPane.h>
#include <Inventor/Win/nodes/SoGuiPosition.h>

// *************************************************************************

SO_NODE_SOURCE(SoGuiPosition);

void
SoGuiPosition::initClass(void)
{
  SO_NODE_INIT_CLASS(SoGuiPosition, SoTransformation, "Transformation");
}

SoGuiPosition::SoGuiPosition(void)
{
  SO_NODE_CONSTRUCTOR(SoGuiPosition);
  SO_NODE_ADD_FIELD(position, (SbVec3f(0.0f, 0.0f, 0.0f)));
}

SoGuiPosition::~SoGuiPosition(void)
{
}

void
SoGuiPosition::doAction(SoAction * action)
{
  int i;
  // SoDebugError::postInfo("SoGuiPosition::doAction", "invoked by %s", action->getTypeId().getName().getString());
  SoGuiPane * pane = NULL;
  const SoFullPath * path = (const SoFullPath *) action->getCurPath();
  for ( i = path->getLength() - 1; (i >= 0) && (pane == NULL); i-- ) {
    SoNode * node = path->getNode(i);
    assert(node);
    if ( node->isOfType(SoGuiPane::getClassTypeId()) ) pane = (SoGuiPane *) node;
  }
  if ( pane == NULL ) {
    SoDebugError::postInfo("SoGuiPosition::GLRender", "SoGuiPosition only works below an SoGuiPane node");
    return;
  }
  pane->moveTo(action->getState(), this->position.getValue());
}

void
SoGuiPosition::GLRender(SoGLRenderAction * action)
{
  this->doAction(action);
}

void
SoGuiPosition::pick(SoPickAction * action)
{
  this->doAction(action);
}

void
SoGuiPosition::rayPick(SoRayPickAction * action)
{
  this->doAction(action);
}

void
SoGuiPosition::getMatrix(SoGetMatrixAction * action)
{
  SoDebugError::postInfo("SoGuiPosition::getMatrix", "invoked");
  int i;
  SoGuiPane * pane = NULL;
  const SoFullPath * path = (const SoFullPath *) action->getCurPath();
  for ( i = path->getLength() - 1; (i >= 0) && (pane == NULL); i-- ) {
    SoNode * node = path->getNode(i);
    assert(node);
    if ( node->isOfType(SoGuiPane::getClassTypeId()) ) pane = (SoGuiPane *) node;
  }
  if ( pane == NULL ) {
    SoDebugError::postInfo("SoGuiPosition::getMatrix", "SoGuiPosition only works below an SoGuiPane node");
    return;
  }
  pane->applyMoveTo(action, this->position.getValue());
}

