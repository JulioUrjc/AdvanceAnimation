// src\Inventor\Win\nodes\Slider1.cpp.  Generated from Slider1.cpp.in by configure.

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <Inventor/errors/SoDebugError.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoHandleEventAction.h>
#include <Inventor/misc/SoState.h>
#include <Inventor/elements/SoLazyElement.h>
#include <Inventor/sensors/SoFieldSensor.h>
#include <Inventor/events/SoMouseButtonEvent.h>
#include <Inventor/events/SoLocation2Event.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/SoLists.h>

#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoTexture2.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoTextureCoordinate2.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoIndexedLineSet.h>

#include <Inventor/Win/common/gl.h>
#include <Inventor/Win/SoWinBasic.h>

#include <Inventor/Win/nodes/SoGuiPane.h>
#include <Inventor/Win/nodes/SoGuiSlider1.h>
#include <assert.h>

// *************************************************************************

/*!
  \class SoGuiSlider1 Inventor/Win/nodes/SoGuiSlider1.h
  \brief A GUI component for a 1-dimensional slider.

  The SoGuiSlider1 node is for creating 2D user interfaces with
  sliders.

  fields:
    min
    max
    value
    orientation - not supported yet
    size
*/

// *************************************************************************

class Slider1 {
public:
  SoGuiSlider1 * kit;
  SoFieldSensor * sizeSensor;
  SoFieldSensor * valueSensor;
  SoFieldSensor * minSensor;
  SoFieldSensor * maxSensor;

  SbBool grabbing;
  float grabpos;
  float grabval;
  float graboffset;
  float pickpos;
  SoGuiPane * pane;

  Slider1(void);
  ~Slider1(void);

  SbColor mincolor, maxcolor;

  // sensors callbacks
  static void sizeChangeCB(void * closure, SoSensor * sensor);
  static void valueChangeCB(void * closure, SoSensor * sensor);
  static void minChangeCB(void * closure, SoSensor * sensor);
  static void maxChangeCB(void * closure, SoSensor * sensor);
};

Slider1::Slider1(void)
{
  this->kit = NULL;
  this->sizeSensor = NULL;
  this->valueSensor = NULL;
  this->minSensor = NULL;
  this->maxSensor = NULL;
  this->grabbing = FALSE;
  this->pane = NULL;
}

#define DELETE_SENSOR(sensor) \
  if ( (sensor) != NULL ) { \
    (sensor)->detach(); \
    delete (sensor); \
    (sensor) = NULL; \
  }

Slider1::~Slider1(void)
{
  DELETE_SENSOR(this->sizeSensor);
  DELETE_SENSOR(this->valueSensor);
  DELETE_SENSOR(this->minSensor);
  DELETE_SENSOR(this->maxSensor);
  this->kit = NULL;
}

#undef DELETE_SENSOR

void
Slider1::sizeChangeCB(void * closure, SoSensor * sensor)
{
  assert(closure);
  Slider1 * internals = (Slider1 *) closure;
  assert(internals->kit);
  internals->kit->sizeUpdate();
}

void
Slider1::valueChangeCB(void * closure, SoSensor * sensor)
{
  assert(closure);
  Slider1 * internals = (Slider1 *) closure;
  assert(internals->kit);
  internals->kit->valueUpdate();
}

void
Slider1::minChangeCB(void * closure, SoSensor * sensor)
{
  assert(closure);
  Slider1 * internals = (Slider1 *) closure;
  assert(internals->kit);
  internals->kit->minUpdate();
}

void
Slider1::maxChangeCB(void * closure, SoSensor * sensor)
{
  assert(closure);
  Slider1 * internals = (Slider1 *) closure;
  assert(internals->kit);
  internals->kit->maxUpdate();
}

// *************************************************************************

#define PRIVATE(obj) ((Slider1 *) obj->internals)

void
SoGuiSlider1::initClass(void)
{
  SO_KIT_INIT_CLASS(SoGuiSlider1, SoBaseKit, "BaseKit");
}

SO_KIT_SOURCE(SoGuiSlider1);

SoGuiSlider1::SoGuiSlider1(void)
{
  this->internals = (void *) new Slider1;
  PRIVATE(this)->kit = this;

  SO_KIT_CONSTRUCTOR(SoGuiSlider1);

  SO_KIT_ADD_FIELD(size, (SbVec3f(1.0f, 1.0f, 0.0f)));
  SO_KIT_ADD_FIELD(orientation, (SoGuiSlider1::X));
  SO_KIT_ADD_FIELD(min, (0.0f));
  SO_KIT_ADD_FIELD(max, (1.0f));
  SO_KIT_ADD_FIELD(value, (0.0f));
  SO_KIT_ADD_FIELD(alwaysHook, (TRUE));

  SO_KIT_DEFINE_ENUM_VALUE(Orientation, X);
  SO_KIT_DEFINE_ENUM_VALUE(Orientation, Y);

  SO_KIT_SET_SF_ENUM_TYPE(orientation, Orientation);

  SO_KIT_ADD_CATALOG_ENTRY(knobLightLineSet, SoIndexedLineSet, FALSE, knobGeometry, "", FALSE);
  SO_KIT_ADD_CATALOG_ENTRY(knobLightMaterial, SoMaterial, FALSE, knobGeometry, knobLightLineSet, TRUE);
  SO_KIT_ADD_CATALOG_ENTRY(knobShadeLineSet, SoIndexedLineSet, FALSE, knobGeometry, knobLightMaterial, FALSE);
  SO_KIT_ADD_CATALOG_ENTRY(knobShadeMaterial, SoMaterial, FALSE, knobGeometry, knobShadeLineSet, TRUE);
  SO_KIT_ADD_CATALOG_ENTRY(knobFaceSet, SoIndexedFaceSet, FALSE, knobGeometry, knobShadeMaterial, FALSE);
  SO_KIT_ADD_CATALOG_ENTRY(knobMaterial, SoMaterial, FALSE, knobGeometry, knobFaceSet, TRUE);
  SO_KIT_ADD_CATALOG_ENTRY(knobCoords, SoCoordinate3, FALSE, knobGeometry, knobMaterial, FALSE);

  SO_KIT_ADD_CATALOG_ENTRY(knobGeometry, SoSeparator, FALSE, topSeparator, "", FALSE);
  SO_KIT_ADD_CATALOG_ENTRY(surfaceFaceSet, SoIndexedFaceSet, FALSE, surfaceGeometry, "", FALSE);
  SO_KIT_ADD_CATALOG_ENTRY(surfaceCoords, SoCoordinate3, FALSE, surfaceGeometry, surfaceFaceSet, FALSE);
  SO_KIT_ADD_CATALOG_ENTRY(surfaceTexCoords, SoTextureCoordinate2, FALSE, surfaceGeometry, surfaceCoords, FALSE);
  SO_KIT_ADD_CATALOG_ENTRY(surfaceTexture, SoTexture2, TRUE, surfaceGeometry, surfaceTexCoords, TRUE);
  SO_KIT_ADD_CATALOG_ENTRY(surfaceMaterial, SoMaterial, TRUE, surfaceGeometry, surfaceTexture, TRUE);
  SO_KIT_ADD_CATALOG_ENTRY(surfaceGeometry, SoSeparator, FALSE, topSeparator, knobGeometry, FALSE);
  SO_KIT_ADD_CATALOG_ENTRY(topSeparator, SoSeparator, FALSE, this, "", FALSE);

  SO_KIT_INIT_INSTANCE();

  static float surfacetexturecoordinates[][2] = { {0.0f, 0.0f}, {1.0f, 0.0f}, {1.0f, 0.0f}, {0.0f, 0.0f} };
  SoTextureCoordinate2 * surfacetexcoords = SO_GET_ANY_PART(this, "surfaceTexCoords", SoTextureCoordinate2);
  assert(surfacetexcoords);
  surfacetexcoords->point.setValues(0, 4, surfacetexturecoordinates);

  static int32_t surfaceindices[] = { 0, 1, 2, -1, 0, 2, 3, -1 };
  SoIndexedFaceSet * surfacefaceset = SO_GET_ANY_PART(this, "surfaceFaceSet", SoIndexedFaceSet);
  assert(surfacefaceset);
  surfacefaceset->textureCoordIndex.setValues(0, 8, surfaceindices);
  surfacefaceset->coordIndex.setValues(0, 8, surfaceindices);
 
  static int32_t knobindices[] = {
    0, 1, 2, -1, 0, 2, 3, -1,
    4, 5, 6, -1, 4, 6, 7, -1,
    8, 2, 5, -1, 8, 5, 9, -1,
    3, 11, 10, -1, 3, 10, 4, -1
  };
  SoIndexedFaceSet * knobfaceset = SO_GET_ANY_PART(this, "knobFaceSet", SoIndexedFaceSet);
  assert(knobfaceset);
  knobfaceset->coordIndex.setValues(0, sizeof(knobindices)/sizeof(knobindices[0]), knobindices);

  this->sizeUpdate();

  SoMaterial * knobmaterial = SO_GET_ANY_PART(this, "knobMaterial", SoMaterial);
  assert(knobmaterial);
  knobmaterial->ambientColor.setValue(0.6f, 0.6f, 0.6f);
  knobmaterial->diffuseColor.setValue(0.6f, 0.6f, 0.6f);
  knobmaterial->emissiveColor.setValue(0.6f, 0.6f, 0.6f);

  SoMaterial * knoblightmaterial = SO_GET_ANY_PART(this, "knobLightMaterial", SoMaterial);
  assert(knoblightmaterial);
  knoblightmaterial->ambientColor.setValue(0.75f, 0.75f, 0.75f);
  knoblightmaterial->diffuseColor.setValue(0.75f, 0.75f, 0.75f);
  knoblightmaterial->emissiveColor.setValue(0.75f, 0.75f, 0.75f);

  SoMaterial * knobshadowmaterial = SO_GET_ANY_PART(this, "knobShadeMaterial", SoMaterial);
  assert(knobshadowmaterial);
  knobshadowmaterial->ambientColor.setValue(0.4f, 0.4f, 0.4f);
  knobshadowmaterial->diffuseColor.setValue(0.4f, 0.4f, 0.4f);
  knobshadowmaterial->emissiveColor.setValue(0.4f, 0.4f, 0.4f);

  // FIXME: move these to correct coordinates
  SoIndexedLineSet * lightlineset = SO_GET_ANY_PART(this, "knobLightLineSet", SoIndexedLineSet);
  assert(lightlineset);
  static int32_t lightindices[] = { 16, 17, 18, -1, 12, 15, 14, -1 };
  lightlineset->coordIndex.setValues(0, sizeof(lightindices) / sizeof(lightindices[0]), lightindices);

  SoIndexedLineSet * shadelineset = SO_GET_ANY_PART(this, "knobShadeLineSet", SoIndexedLineSet);
  assert(shadelineset);
  static int32_t shadeindices[] = { 12, 13, 14, -1, 16, 19, 18, -1 };
  shadelineset->coordIndex.setValues(0, sizeof(shadeindices) / sizeof(shadeindices[0]), shadeindices);

  // set up sensors
  PRIVATE(this)->sizeSensor = new SoFieldSensor(Slider1::sizeChangeCB, PRIVATE(this));
  PRIVATE(this)->sizeSensor->attach(&(this->size));
  PRIVATE(this)->valueSensor = new SoFieldSensor(Slider1::valueChangeCB, PRIVATE(this));
  PRIVATE(this)->valueSensor->attach(&(this->value));
  PRIVATE(this)->minSensor = new SoFieldSensor(Slider1::minChangeCB, PRIVATE(this));
  PRIVATE(this)->minSensor->attach(&(this->min));
  PRIVATE(this)->maxSensor = new SoFieldSensor(Slider1::maxChangeCB, PRIVATE(this));
  PRIVATE(this)->maxSensor->attach(&(this->max));
}

SoGuiSlider1::~SoGuiSlider1(void)
{
  Slider1 * obj = PRIVATE(this);
  delete obj;
  this->internals = NULL;
}

void
SoGuiSlider1::setSurfaceColor(const SbColor & valuearg)
{
  // FIXME: use Material or basecolor instead of texture
  this->setSurfaceColor(valuearg, valuearg);
#if 0
  PRIVATE(this)->mincolor = valuearg;
  PRIVATE(this)->maxcolor = valuearg;
  this->setPart("surfaceTexture", NULL);
#endif
}

void
SoGuiSlider1::setSurfaceColor(const SbColor & minvalue, const SbColor & maxvalue)
{
  PRIVATE(this)->mincolor = minvalue;
  PRIVATE(this)->maxcolor = maxvalue;

  this->setPart("surfaceMaterial", NULL);
  SoTexture2 * texturenode = SO_GET_ANY_PART(this, "surfaceTexture", SoTexture2);
  assert(texturenode);

  texturenode->image.setValue(SbVec2s(256, 1), 3, NULL);
  texturenode->model.setValue(SoTexture2::DECAL);

  SbVec2s sizeval;
  int nc;
  unsigned char * buf = texturenode->image.startEditing(sizeval, nc);
  float rmin = minvalue[0];
  float gmin = minvalue[1];
  float bmin = minvalue[2];
  float rmax = maxvalue[0];
  float gmax = maxvalue[1];
  float bmax = maxvalue[2];
  for ( int x = 0; x < sizeval[0]; x += 1 ) {
    buf[x*nc+0] = (unsigned char) ((rmin + ((float) x / (float) (sizeval[0]-1)) * (rmax - rmin)) * 255.0f);
    buf[x*nc+1] = (unsigned char) ((gmin + ((float) x / (float) (sizeval[0]-1)) * (gmax - gmin)) * 255.0f);
    buf[x*nc+2] = (unsigned char) ((bmin + ((float) x / (float) (sizeval[0]-1)) * (bmax - bmin)) * 255.0f);
    for ( int y = 1; y < sizeval[1]; y += 1 ) {
      buf[(y*sizeval[0]+x)*nc+0] = buf[x*nc+0];
      buf[(y*sizeval[0]+x)*nc+1] = buf[x*nc+1];
      buf[(y*sizeval[0]+x)*nc+2] = buf[x*nc+2];
    }
  }
  texturenode->image.finishEditing();
}

SbColor
SoGuiSlider1::getValueAsColor(void) const
{
  // FIXME: support custom textures
  float val = this->value.getValue();
  float minval = this->min.getValue();
  float maxval = this->max.getValue();
  float factor = (maxval - minval) / (val - minval);
  float r = SoWinClamp(PRIVATE(this)->mincolor[0] + (PRIVATE(this)->maxcolor[0] - PRIVATE(this)->mincolor[0]) * factor, 0.0f, 1.0f);
  float g = SoWinClamp(PRIVATE(this)->mincolor[1] + (PRIVATE(this)->maxcolor[1] - PRIVATE(this)->mincolor[1]) * factor, 0.0f, 1.0f);
  float b = SoWinClamp(PRIVATE(this)->mincolor[2] + (PRIVATE(this)->maxcolor[2] - PRIVATE(this)->mincolor[2]) * factor, 0.0f, 1.0f);
  return SbColor(r, g, b);
}

void
SoGuiSlider1::sizeUpdate(void)
{
  SbVec3f sizeval = this->size.getValue();
  if ( sizeval[0] != 0.0f && sizeval[1] != 0.0f ) {
    float coordinates[][3] = { {0.0f, 0.0f, 0.0f}, {sizeval[0], 0.0f, 0.0f}, {sizeval[0], sizeval[1], 0.0f}, {0.0f, sizeval[1], 0.0f} };
    SoCoordinate3 * coords = SO_GET_ANY_PART(this, "surfaceCoords", SoCoordinate3);
    assert(coords);
    coords->point.setValues(0, sizeof(coordinates) / sizeof(coordinates[0]), coordinates);
    this->valueUpdate();
  }
}

void
SoGuiSlider1::valueUpdate(void)
{
  SbVec3f sizeval = this->size.getValue();
  float val = this->value.getValue();
  float minval = this->min.getValue();
  float maxval = this->max.getValue();
  if ( minval < maxval ) {
    if ( val < minval ) {
      this->value.setValue(minval);
      val = minval;
    } else if ( val > maxval ) {
      this->value.setValue(maxval);
      val = maxval;
    }
  } else {
    // we also support inverse sliders where min > max
    if ( val > minval ) {
      this->value.setValue(minval);
      val = minval;
    } else if ( val < maxval ) {
      this->value.setValue(maxval);
      val = maxval;
    }
  }
  // store previous height & value to avoid redundant updates
  float voff = (float) floor(((val - minval) / (maxval - minval)) * sizeval[0]);
  float knobcoordinates[][3] = {
    // faces
    {-7.0f+voff, -4.0f, 0.0f}, {8.0f+voff, -4.0f, 0.0f}, {8.0f+voff, -1.0f, 0.0f}, {-7.0f+voff, -1.0f, 0.0f}, 
    {-7.0f+voff, sizeval[1]+1.0f, 0.0f}, {8.0f+voff, sizeval[1]+1.0f, 0.0f}, {8.0f+voff, sizeval[1]+4.0f, 0.0f}, {-7.0f+voff, sizeval[1]+4.0f, 0.0f}, 
    {3.0f+voff, -1.0f, 0.0f}, {3.0f+voff, sizeval[1]+1.0f, 0.0f}, {-2.0f+voff, sizeval[1]+1.0f, 0.0f}, {-2.0f+voff, -1.0f, 0.0f},
    // outside lines
    // also used in handleEvent() so don't change their significance based on index
    {-8.0f+voff, -5.0f, 0.0f}, {8.0f+voff, -5.0f, 0.0f}, {8.0f+voff, sizeval[1]+4.0f, 0.0f}, {-8.0f+voff, sizeval[1]+4.0f, 0.0f},
    // inside lines
    {-2.0f+voff, -1.0f, 0.0f}, {2.0f+voff, -1.0f, 0.0f}, {2.0f+voff, sizeval[1], 0.0f}, {-2.0f+voff, sizeval[1], 0.0f}
  };
  SoCoordinate3 * knobcoords = SO_GET_ANY_PART(this, "knobCoords", SoCoordinate3);
  assert(knobcoords);
  knobcoords->point.setValues(0, sizeof(knobcoordinates) / sizeof(knobcoordinates[0]), knobcoordinates);
}

void
SoGuiSlider1::minUpdate(void)
{
  float minval = this->min.getValue();
  float maxval = this->max.getValue();
  float val = this->value.getValue();
  if ( minval < maxval ) {
    if ( val < minval ) this->value.setValue(minval);
  } else {
    if ( val > minval ) this->value.setValue(minval);
  }
}

void
SoGuiSlider1::maxUpdate(void)
{
  float minval = this->min.getValue();
  float maxval = this->max.getValue();
  float val = this->value.getValue();
  if ( minval < maxval ) {
    if ( val > maxval ) this->value.setValue(maxval);
  } else {
    if ( val < maxval ) this->value.setValue(maxval);
  }
}

void
SoGuiSlider1::handleEvent(SoHandleEventAction * action)
{
  if ( action->isHandled() ) return;
  const SoEvent * event = action->getEvent();

  if ( PRIVATE(this)->grabbing ) { // click-and-drag
    if ( event->isOfType(SoLocation2Event::getClassTypeId()) ) {
      assert(PRIVATE(this)->pane != NULL);
      // although the return value is discarded, we need to make this call to make
      // sure the raypick action is run over the scene graph so the pane can return
      // a useful value...
      action->getPickedPoint();
      SbVec2f raypos = PRIVATE(this)->pane->getRayPickIntersectionPoint();
      if ( raypos[0] != -1.0f ) {
        float imagpos = raypos[0] + PRIVATE(this)->graboffset;
        float minval = this->min.getValue();
        float maxval = this->max.getValue();
        float imagval;
        if ( minval < maxval ) imagval = SoWinClamp(imagpos / this->size.getValue()[0], minval, maxval);
        else imagval = SoWinClamp(imagpos / this->size.getValue()[0], maxval, minval);
        this->value.setValue(imagval);
      }
      action->setHandled();
    }
    else if ( event->isOfType(SoMouseButtonEvent::getClassTypeId()) ) {
      SoMouseButtonEvent * mbevent = (SoMouseButtonEvent *) event;
      if ( (mbevent->getButton() == SoMouseButtonEvent::BUTTON1) &&
           (mbevent->getState() == SoButtonEvent::UP) ) {
        PRIVATE(this)->grabbing = FALSE;
        PRIVATE(this)->pane = NULL;
        action->setHandled();
      }
    }
  }
  else if ( event->isOfType(SoMouseButtonEvent::getClassTypeId()) ) {
    SoMouseButtonEvent * mbevent = (SoMouseButtonEvent *) event;
    if ( (mbevent->getButton() == SoMouseButtonEvent::BUTTON1) &&
         (mbevent->getState() == SoButtonEvent::DOWN) ) {
      action->setPickRadius(0);
      const SoPickedPointList & pplist = action->getPickedPointList();
      if ( pplist.getLength() > 0 ) {
        int i;
        for ( i = 0; i < pplist.getLength(); i++ ) {
          if ( action->isHandled() ) break;
          const SoPickedPoint * pp = pplist[i];
          const SoPath * path = pp->getPath();
          SoNode * node = ((SoFullPath *) path)->getTail();
          if ( node == ((SoNode *) SO_GET_ANY_PART(this, "knobFaceSet", SoIndexedFaceSet)) ) {
            SbVec3f point = pp->getObjectPoint();
            SbVec3f sizeval = this->size.getValue();
            SoCoordinate3 * knobcoords = SO_GET_ANY_PART(this, "knobCoords", SoCoordinate3);
            assert(knobcoords);
            SbVec3f knobmin = knobcoords->point[12];
            SbVec3f knobmax = knobcoords->point[14];
            if ( point[0] >= knobmin[0] && point[0] <= knobmax[0] &&
                 point[1] >= knobmin[1] && point[1] <= knobmax[1] ) {
              PRIVATE(this)->grabbing = TRUE;
              action->setHandled();
  
              const SoFullPath * path = (const SoFullPath *) action->getCurPath();
              int i = path->getLength() - 1;
              SoNode * node = NULL;
              for ( ; i >= 0; i-- ) {
                node = path->getNode(i);
                if ( node->isOfType(SoGuiPane::getClassTypeId()) ) break;
                node = NULL;
              }
              assert(node != NULL);
              PRIVATE(this)->pane = (SoGuiPane *) node;
              PRIVATE(this)->grabval = this->value.getValue();
              PRIVATE(this)->pickpos = point[0];
  
              SbVec2f raypos = PRIVATE(this)->pane->getRayPickIntersectionPoint();
              PRIVATE(this)->grabpos = raypos[0];
              float realval = ((this->value.getValue() - this->min.getValue()) /
                (this->max.getValue() - this->min.getValue())) * this->size.getValue()[0];
              PRIVATE(this)->graboffset = realval - raypos[0];
            }
          }
        }
        for ( i = 0; i < pplist.getLength(); i++ ) {
          if ( action->isHandled() ) break;
          const SoPickedPoint * pp = pplist[i];
          const SoPath * path = pp->getPath();
          SoNode * node = ((SoFullPath *) path)->getTail();
          if ( node == ((SoNode *) SO_GET_ANY_PART(this, "surfaceFaceSet", SoIndexedFaceSet)) ) {
            SbVec3f point = pp->getObjectPoint();
            SbVec3f sizeval = this->size.getValue();
            this->value = this->min.getValue() + ((point[0] / sizeval[0]) * (this->max.getValue() - this->min.getValue()));
            action->setHandled();
            if ( this->alwaysHook.getValue() ) {
              PRIVATE(this)->grabbing = TRUE;
              const SoFullPath * path = (const SoFullPath *) action->getCurPath();
              int i = path->getLength() - 1;
              SoNode * node = NULL;
              for ( ; i >= 0; i-- ) {
                node = path->getNode(i);
                if ( node->isOfType(SoGuiPane::getClassTypeId()) ) break;
                node = NULL;
              }
              assert(node != NULL);
              PRIVATE(this)->pane = (SoGuiPane *) node;
              PRIVATE(this)->grabval = this->value.getValue();
              PRIVATE(this)->pickpos = point[0];
              SbVec2f raypos = PRIVATE(this)->pane->getRayPickIntersectionPoint();
                PRIVATE(this)->grabpos = raypos[0];
              float realval = ((this->value.getValue() - this->min.getValue()) /
                (this->max.getValue() - this->min.getValue())) * this->size.getValue()[0];
              PRIVATE(this)->graboffset = realval - raypos[0];
            }
          }
        }
      }
    }
  }
}

#undef PRIVATE

