/* ================================================================
 *
 * PyCA Project
 *
 * Copyright (c) J. Samuel Preston, Linh K. Ha, Sarang C. Joshi. All
 * rights reserved.  See Copyright.txt or for details.
 *
 * This software is distributed WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the above copyright notice for more information.
 *
 * ================================================================ */

#ifndef __SELECTOR__H
#define __SELECTOR__H

namespace PyCA {

template<int id, typename U, typename U1, typename U2>
struct Selector{
    typedef U Result;
};

template<typename U, typename U1, typename U2>
struct Selector<1, U, U1, U2>{
    typedef U1 Result;
};

template<typename U, typename U1, typename U2>
struct Selector<2, U, U1, U2>{
    typedef U2 Result;
};


template<int id, typename U, typename U1>
struct BinSelector{
    typedef U Result;
};

template<typename U, typename U1>
struct BinSelector<1, U, U1>{
    typedef U1 Result;
};

} // end namespace PyCA

#endif
