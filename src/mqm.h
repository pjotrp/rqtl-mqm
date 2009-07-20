/**********************************************************************
 *
 * mqm.h
 *
 * Master include file for common mqm.headers
 *
 * copyright (c) 2009 Ritsert Jansen, Danny Arends, Pjotr Prins and Karl Broman
 *
 * last modified July, 2009
 * first written July, 2009
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License,
 *     version 3, as published by the Free Software Foundation.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the GNU
 *     General Public License, version 3, for more details.
 *
 *     A copy of the GNU General Public License, version 3, is available
 *     at http://www.r-project.org/Licenses/GPL-3
 *
 **********************************************************************/


#ifndef __MQM_H
  #define __MQM_H

  #include <R.h>

  // #include <R_ext/PrtUtil.h>
  // #include <R_ext/RS.h> 
  // #include <R_ext/Utils.h>
  #include "standalone.h"
  #include "util.h"
  #include "mqmdatatypes.h"
  #include "mqmprob.h"        
  #include "mqmmixture.h"    
  #include "mqmregression.h"
  #include "mqmaugment.h"
  #include "mqmeliminate.h"
  #include "mqmmapqtl.h"  
  #include "mqmscan.h"

const unsigned char MLEFT     ='L';
const unsigned char MRIGHT    ='R';
const unsigned char MMIDDLE   ='M';
const unsigned char MUNLINKED ='U';
 
const unsigned char MA       = '0';  // Homozygous parent A
const unsigned char MB       = '1';  // Homyzygous parent B
const unsigned char MH       = '2';  // Heterozygous
const unsigned char MNOTA    = '3';  // Not A
const unsigned char MNOTB    = '4';  // Not B 
const unsigned char MMISSING = '9';  // Uknown

#endif
