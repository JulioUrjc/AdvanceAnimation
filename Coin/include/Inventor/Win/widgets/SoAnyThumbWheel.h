#ifndef SOANY_THUMBWHEEL_H
#define SOANY_THUMBWHEEL_H

// src\Inventor\Win\widgets\SoAnyThumbWheel.h.  Generated from SoAnyThumbWheel.h.in by configure.

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

class SoAnyThumbWheel {
public:
  enum State              { DISABLED, ENABLED };
  enum Alignment          { VERTICAL, HORIZONTAL };
  enum BoundaryHandling   { MODULATE, ACCUMULATE, CLAMP };
  enum Movement           { UNIFORM, AUTHENTIC };
  enum GraphicsByteOrder  { ABGR, RGBA, ARGB, BGRA };

  SoAnyThumbWheel(void);
  ~SoAnyThumbWheel(void);

  void setSize(const int diameter, const int width);
  void getSize(int & diameter, int & width) const;

  void setColor(const float red, const float green, const float blue);
  void getColor(float & red, float & green, float & blue) const;
  void setColorFactors(const float light, const float front, const float normal, const float shade);
  void getColorFactors(float & light, float & front, float & normal, float & shade) const;

  int getNumBitmaps(void) const;
  void drawBitmap(const int number, void * bitmap, Alignment alignment) const;
  float calculateValue(const float origValue, const int origPosition, const int deltaPosition) const;
  int getBitmapForValue(const float value, const State state) const;

  void setGraphicsByteOrder(const GraphicsByteOrder byteorder);
  GraphicsByteOrder getGraphicsByteOrder(void) const;

  void setMovement(const Movement movement);
  Movement getMovement(void) const;

  void setBoundaryHandling(const BoundaryHandling handling);
  BoundaryHandling getBoundaryHandling(void) const;

private:
  unsigned int swapWord(unsigned int) const;

  int diameter, width;
  // float disabledred, disabledgreen, disabledblue; // not implemented
  float red, green, blue;
  float light, front, normal, shade;

  GraphicsByteOrder  byteorder;
  BoundaryHandling   boundaryhandling;
  Movement           movement;

  enum Tables { SIN, COS, RAD, NUMTABLES };

  mutable float * tables [ NUMTABLES ];
  mutable int dirtyTables;
  mutable int dirtyVariables;
  mutable float squarelength, squarespacing, shadelength, unistep, numsquares;

  void drawDisabledWheel(const int number, void * bitmap, Alignment alignment) const;
  void drawEnabledWheel(const int number, void * bitmap, Alignment alignment) const;

  void validate(void) const;

}; // class SoAnyThumbWheel

// ************************************************************************

#endif // ! SOANY_THUMBWHEEL_H
