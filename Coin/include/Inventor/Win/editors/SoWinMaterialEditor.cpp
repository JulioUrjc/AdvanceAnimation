// src\Inventor\Win\editors\SoWinMaterialEditor.cpp.  Generated from SoGuiMaterialEditor.cpp.in by configure.

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
#endif // HAVE_CONFIG_H

#include <Inventor/SbPList.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/sensors/SoNodeSensor.h>
#ifdef HAVE_SOVRMLMATERIAL
#include <Inventor/VRMLnodes/SoVRMLMaterial.h>
#endif // HAVE_SOVRMLMATERIAL

#include <Inventor/Win/SoAny.h>
#include <Inventor/Win/nodes/SoGuiMaterialEditor.h>
#include <Inventor/Win/editors/SoWinColorEditor.h>
#include <Inventor/Win/editors/SoWinMaterialEditor.h>

// *************************************************************************

class MaterialEditorComponent {
public:
  MaterialEditorComponent(void);
  ~MaterialEditorComponent(void);

  SoWinMaterialEditor * api;

  void constructor(void);

  SoWinMaterialEditor::UpdateFrequency frequency;

  SbPList * callbacks;
  void invokeCallbacks(void);

  SoNodeSensor * attachment;
  SoMaterial * material;
  SoGuiMaterialEditor * editor;

  static SoGuiColorEditor * coloreditoropen_cb(void * closure);
  static void coloreditorclose_cb(void * closure, SoGuiColorEditor * editor);
  SoWinColorEditor * coloreditor;

  static SoSeparator * getSuperScene(void);
  static SoNode * getSceneGraph(void);
  static const char * superscene[];
};

// *************************************************************************

#define PRIVATE(obj) ((MaterialEditorComponent *) obj->internals)

SOWIN_OBJECT_SOURCE(SoWinMaterialEditor);

SoWinMaterialEditor::SoWinMaterialEditor(HWND parent, const char * name, SbBool embed)
  : inherited(parent, name, embed)
{
  this->internals = new MaterialEditorComponent;
  PRIVATE(this)->api = this;
  PRIVATE(this)->constructor();
}

SoWinMaterialEditor::SoWinMaterialEditor(HWND parent, const char * const name, SbBool embed, SbBool build)
: inherited(parent, name, embed)
{
  this->internals = new MaterialEditorComponent;
  PRIVATE(this)->api = this;
  PRIVATE(this)->constructor();
}

SoWinMaterialEditor::~SoWinMaterialEditor(void)
{
  assert(PRIVATE(this)->editor);
  MaterialEditorComponent * obj = PRIVATE(this);
  delete obj;
}

void
SoWinMaterialEditor::attach(SoMaterial * material, int index)
{
  assert(PRIVATE(this)->editor);
  PRIVATE(this)->editor->material.setValue(material);
  PRIVATE(this)->editor->index.setValue(index);
}

void
SoWinMaterialEditor::attach(SoVRMLMaterial * material)
{
  assert(PRIVATE(this)->editor);
#ifdef HAVE_SOVRMLMATERIAL
  PRIVATE(this)->editor->material.setValue(material);
#else // !HAVE_SOVRMLMATERIAL
  PRIVATE(this)->editor->material.setValue((SoNode *) material);
#endif // !HAVE_SOVRMLMATERIAL
}

void
SoWinMaterialEditor::detach(void)
{
  PRIVATE(this)->editor->material.setValue(NULL);
}

SbBool
SoWinMaterialEditor::isAttached(void)
{
  assert(PRIVATE(this)->editor != NULL);
  return (PRIVATE(this)->editor->material.getValue() != NULL) ? TRUE : FALSE;
}

void
SoWinMaterialEditor::addMaterialChangedCallback(SoWinMaterialEditorCB * callback, void * closure)
{
  if ( PRIVATE(this)->callbacks == NULL )
    PRIVATE(this)->callbacks = new SbPList;
  PRIVATE(this)->callbacks->append((void *) callback);
  PRIVATE(this)->callbacks->append(closure);
}

void
SoWinMaterialEditor::removeMaterialChangedCallback(SoWinMaterialEditorCB * callback, void * closure)
{
  assert(PRIVATE(this)->callbacks);
  int i;
  for ( i = PRIVATE(this)->callbacks->getLength(); i > 0; i -= 2 ) {
    SoWinMaterialEditorCB * cb = (SoWinMaterialEditorCB *) (*(PRIVATE(this)->callbacks))[i-2];
    void * cl = (*(PRIVATE(this)->callbacks))[i-1];
    if ( cb == callback && cl == closure ) {
      PRIVATE(this)->callbacks->remove(i-1);
      PRIVATE(this)->callbacks->remove(i-2);
    }
  }
}

void
SoWinMaterialEditor::setUpdateFrequency(SoWinMaterialEditor::UpdateFrequency frequency)
{
  PRIVATE(this)->frequency = frequency;
}

SoWinMaterialEditor::UpdateFrequency
SoWinMaterialEditor::getUpdateFrequency(void) const
{
  return PRIVATE(this)->frequency;
}

void
SoWinMaterialEditor::setMaterial(const SoMaterial & material)
{
  PRIVATE(this)->material->ambientColor.setValue(material.ambientColor[0]);
  PRIVATE(this)->material->diffuseColor.setValue(material.diffuseColor[0]);
  PRIVATE(this)->material->specularColor.setValue(material.specularColor[0]);
  PRIVATE(this)->material->emissiveColor.setValue(material.emissiveColor[0]);
  PRIVATE(this)->material->shininess.setValue(material.shininess[0]);
  PRIVATE(this)->material->transparency.setValue(material.transparency[0]);
}

void
SoWinMaterialEditor::setMaterial(const SoVRMLMaterial & material)
{
#ifdef HAVE_SOVRMLMATERIAL
  float h, s, v;
  h = 0.0f;
  s = 0.0f;
  v = material.ambientIntensity.getValue();
  SbColor ambient;
  ambient.setHSVValue(h, s, v);
  PRIVATE(this)->material->ambientColor.setValue(ambient);
  PRIVATE(this)->material->diffuseColor.setValue(material.diffuseColor.getValue());
  PRIVATE(this)->material->specularColor.setValue(material.specularColor.getValue());
  PRIVATE(this)->material->emissiveColor.setValue(material.emissiveColor.getValue());
  PRIVATE(this)->material->shininess.setValue(material.shininess.getValue());
  PRIVATE(this)->material->transparency.setValue(material.transparency.getValue());
#endif // HAVE_SOVRMLMATERIAL
}

const SoMaterial &
SoWinMaterialEditor::getMaterial(void) const
{
  return *(PRIVATE(this)->material);
}

SbBool
SoWinMaterialEditor::isAttachedVRML(void)
{
  SoNode * material = PRIVATE(this)->editor->material.getValue();
#ifdef HAVE_SOVRMLMATERIAL
  return (material && material->isOfType(SoVRMLMaterial::getClassTypeId())) ? TRUE : FALSE;
#endif
  return FALSE;
}

SoGuiMaterialEditor *
SoWinMaterialEditor::getEditor(void) const
{
  assert(0);
  return NULL;
}

// *************************************************************************

const char *
SoWinMaterialEditor::getDefaultWidgetName(void) const
{
  static const char widgetName[] = "SoWinMaterialEditor";
  return widgetName;
}

const char *
SoWinMaterialEditor::getDefaultTitle(void) const
{
  static const char title[] = "MaterialEditor";
  return title;
}

const char *
SoWinMaterialEditor::getDefaultIconTitle(void) const
{
  static const char iconTitle[] = "MatEd";
  return iconTitle;
}

#undef PRIVATE

// *************************************************************************
// MaterialEditorComponent
// *************************************************************************

#define PUBLIC(obj)  (obj->api)

const char *
MaterialEditorComponent::superscene[] =
{
  "#Inventor V2.1 ascii",
  "",
  "Separator {",
  "  DirectionalLight { direction 0 0 -1 color 1 1 1 intensity 0.8 }",
  "  OrthographicCamera { }",
  "  DEF viewportfix SoGuiViewportFix { }",
  "  Material { ambientColor 0.8 0.8 0.8 }",
  "}",
  NULL
};

// *************************************************************************

MaterialEditorComponent::MaterialEditorComponent(void)
{
  this->api = NULL;
  this->attachment = NULL;
  this->editor = NULL;
  this->callbacks = NULL;
  this->material = NULL;
}

MaterialEditorComponent::~MaterialEditorComponent(void)
{
  if ( this->attachment != NULL ) {
    delete this->attachment;
  }
  if ( this->callbacks != NULL )
    delete this->callbacks;
}

void
MaterialEditorComponent::constructor(void)
{
  PUBLIC(this)->setSize(SbVec2s(550,300));
  SoNode * scene = MaterialEditorComponent::getSceneGraph();
  PUBLIC(this)->setSceneGraph(scene);
  assert(scene->isOfType(SoSeparator::getClassTypeId()));
  SoSeparator * sceneroot = (SoSeparator *) scene;
  this->editor = (SoGuiMaterialEditor *) sceneroot->getChild(sceneroot->getNumChildren()-1);
  assert(this->editor->isOfType(SoGuiMaterialEditor::getClassTypeId()));
  this->material = this->editor->getSphereMaterialNode();
  this->editor->setColorEditorCallbacks(MaterialEditorComponent::coloreditoropen_cb, MaterialEditorComponent::coloreditorclose_cb, this);
}

void
MaterialEditorComponent::invokeCallbacks(void)
{
  // FIXME
}

SoSeparator *
MaterialEditorComponent::getSuperScene(void)
{
  SoNode * root = SoAny::loadSceneGraph(MaterialEditorComponent::superscene);
  assert(root != NULL);
  assert(root->isOfType(SoSeparator::getClassTypeId()));
  return (SoSeparator *) root;
}

SoNode *
MaterialEditorComponent::getSceneGraph(void)
{
  SoSeparator * root = MaterialEditorComponent::getSuperScene();
  assert(root != NULL);
  // FIXME: route ViewportFix field value into editor
  root->addChild(new SoGuiMaterialEditor);
  return root;
}

SoGuiColorEditor *
MaterialEditorComponent::coloreditoropen_cb(void * closure)
{
  MaterialEditorComponent * me = (MaterialEditorComponent *) closure;
  me->coloreditor = new SoWinColorEditor;
  me->coloreditor->show();
  return me->coloreditor->getEditor();
}

void
MaterialEditorComponent::coloreditorclose_cb(void * closure, SoGuiColorEditor * editor)
{
  assert(closure);
  MaterialEditorComponent * me = (MaterialEditorComponent *) closure;
  me->coloreditor->hide();
  delete me->coloreditor;
  me->coloreditor = NULL;
}

#undef PUBLIC

// *************************************************************************
