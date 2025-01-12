// src\Inventor\Win\engines\Format.cpp.  Generated from Format.cpp.in by configure.

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

#ifndef SOWIN_INTERNAL
#error this is a private header file
#endif // !SOWIN_INTERNAL

#include <Inventor/SbString.h>
#include <Inventor/errors/SoDebugError.h>
#include <stdio.h>

#include <Inventor/Win/engines/SoGuiFormat.h>

void
SoGuiFormat::initClass(void)
{
  SO_ENGINE_INIT_CLASS(SoGuiFormat, SoEngine, "Engine");
}

SO_ENGINE_SOURCE(SoGuiFormat);

SoGuiFormat::SoGuiFormat(void)
{
  this->internals = NULL;

  SO_ENGINE_CONSTRUCTOR(SoGuiFormat);

  SO_ENGINE_ADD_INPUT(float1, (0.0f));
  SO_ENGINE_ADD_INPUT(format, (""));

  SO_ENGINE_ADD_OUTPUT(output, SoSFString);
}

SoGuiFormat::~SoGuiFormat(void)
{
}

void
SoGuiFormat::evaluate(void)
{
#ifdef __COIN__
  SbString result;
  result.sprintf(this->format.getValue().getString(), this->float1.getValue());
#else // __COIN__
  char buf[4096]; // FIXME: temporary workaround
  sprintf(buf, this->format.getValue().getString(), this->float1.getValue());
  SbString result = buf;
#endif // ! __COIN__
 
  // SoDebugError::postInfo("SoGuiFormat::evaluate", "%s", result.getString());
  SO_ENGINE_OUTPUT(output, SoSFString, setValue(result));
}

