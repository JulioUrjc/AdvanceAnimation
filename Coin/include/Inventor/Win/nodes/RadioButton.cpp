// src\Inventor\Win\nodes\RadioButton.cpp.  Generated from RadioButton.cpp.in by configure.

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

#include <Inventor/actions/SoHandleEventAction.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/sensors/SoFieldSensor.h>
#include <Inventor/events/SoMouseButtonEvent.h>
#include <Inventor/SoPickedPoint.h>

#include <Inventor/Win/SoAny.h>
#include <Inventor/Win/nodes/SoGuiRadioButton.h>
#include <assert.h>

// *************************************************************************

class RadioButton {
public:
  RadioButton(void);

  SoGuiRadioButton * api;

  SoCoordinate3 * coords;
  SoNode * faceset;
  SoFieldSensor * size_sensor;
  static void size_updated_cb(void * closure, SoSensor * sensor);

  static const char * scene[];
};

// *************************************************************************

#define PRIVATE(obj) ((RadioButton *)obj->internals)

void
SoGuiRadioButton::initClass(void)
{
  SO_KIT_INIT_CLASS(SoGuiRadioButton, SoBaseKit, "BaseKit");
}

SO_KIT_SOURCE(SoGuiRadioButton);

SoGuiRadioButton::SoGuiRadioButton(void)
{
  this->internals = new RadioButton;
  PRIVATE(this)->api = this;

  SO_KIT_CONSTRUCTOR(SoGuiRadioButton);

  SO_KIT_ADD_FIELD(size, (SbVec3f(1.0f, 1.0f, 0.0f)));
  SO_KIT_ADD_FIELD(on, (FALSE));

  SO_KIT_ADD_CATALOG_ENTRY(root, SoSeparator, FALSE, this, "", FALSE);

  SO_KIT_INIT_INSTANCE();

  SoNode * scene = SoAny::loadSceneGraph(RadioButton::scene);
  assert(scene);
  assert(scene->isOfType(SoSeparator::getClassTypeId()));
  scene->ref();

  PRIVATE(this)->coords = (SoCoordinate3 *) SoAny::scanSceneForName(scene, "coords");
  assert(PRIVATE(this)->coords);
  assert(PRIVATE(this)->coords->isOfType(SoCoordinate3::getClassTypeId()));
  PRIVATE(this)->faceset = SoAny::scanSceneForName(scene, "faceset");
  assert(PRIVATE(this)->faceset);

  scene->unrefNoDelete();
  this->setAnyPart("root", scene);

  PRIVATE(this)->size_sensor = new SoFieldSensor(RadioButton::size_updated_cb, PRIVATE(this));
  PRIVATE(this)->size_sensor->attach(&(this->size));
}


SoGuiRadioButton::~SoGuiRadioButton(void)
{
  delete PRIVATE(this)->size_sensor;
  RadioButton * obj = PRIVATE(this);
  delete obj;
}

void
SoGuiRadioButton::handleEvent(SoHandleEventAction * action)
{
  const SoEvent * ev = action->getEvent();
  if ( ev->isOfType(SoMouseButtonEvent::getClassTypeId()) ) {
    SbBool hit = FALSE;
    const SoPickedPointList & ppoints = action->getPickedPointList();
    assert(PRIVATE(this)->faceset);
    int i = 0;
    for ( i = 0; !hit && i < ppoints.getLength(); i++ ) {
      const SoPickedPoint * point = ppoints[i];
      const SoFullPath * path = (const SoFullPath *) point->getPath();
      assert(path);
      SoNode * node = path->getTail();
      if ( node == PRIVATE(this)->faceset ) hit = TRUE;
    }
    if ( hit ) {
      const SoMouseButtonEvent * event = (SoMouseButtonEvent *) ev;
      if ( event->getState() == SoButtonEvent::DOWN ) {
        this->on.setValue(TRUE);
        action->setHandled();
      }
    }
  }
  if ( !action->isHandled() ) {
    inherited::handleEvent(action);
  }
}

#undef PRIVATE

// *************************************************************************
// RadioButton
// *************************************************************************

#define PUBLIC(obj) (obj->api)

const char *
RadioButton::scene[] =
{
  "#Inventor V2.1 ascii",
  "",
  "Separator {",
  "  DEF coords Coordinate3 {",
  "    point [",
  "      0 0 0,",
  "      1 0 0,",
  "      1 1 0,",
  "      0 1 0",
  "    ]",
  "  }",
  "  DEF faceset IndexedFaceSet {",
  "    coordIndex [",
  "      0 1 2 -1",
  "      0 2 3 -1",
  "    ]",
  "  }",
  "}",
  NULL
};

RadioButton::RadioButton(void)
{
  this->api = NULL;
  this->coords = NULL;
  this->size_sensor = NULL;
}

void
RadioButton::size_updated_cb(void * closure, SoSensor * sensor)
{
  assert(closure);
  RadioButton * me = (RadioButton *) closure;
  assert(PUBLIC(me));
  SbVec3f size = PUBLIC(me)->size.getValue();
  assert(me->size_sensor);
  me->size_sensor->detach();
  assert(me->coords);
  SbBool save = me->coords->point.enableNotify(FALSE);
  me->coords->point.set1Value(0, SbVec3f(0.0f, 0.0f, 0.0f));
  me->coords->point.set1Value(1, SbVec3f(size[0], 0.0f, 0.0f));
  me->coords->point.set1Value(2, SbVec3f(size[0], size[1], 0.0f));
  me->coords->point.set1Value(3, SbVec3f(0.0f, size[1], 0.0f));
  me->coords->enableNotify(save);
  if ( save ) me->coords->point.touch();
  me->size_sensor->attach(&(PUBLIC(me)->size));
}

#undef PUBLIC

// *************************************************************************
