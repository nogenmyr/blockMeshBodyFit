/**
 * \file includes/Log.h
 * \brief A logging facility.
 * 
 * Copyright 2009 Tobias Preclik
 * 
 * This file is part of amgpp.
 * 
 * Amgpp is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Amgpp is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with amgpp.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef LOG_H_
#define LOG_H_

#define INFO     0
#define WARNING  1
#define CRITICAL 2
#define ERROR    3

extern unsigned int loglevel_;

#endif
