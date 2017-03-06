/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
/**
 * @file Granulometry.cpp
 * @ingroup Tools
 * @author Nicolas Normand (\c Nicolas.Normand@polytech.univ-nantes.fr)
 * Université Bretagne Loire, Université de Nantes,
 * Laboratoire des Sciences du Numérique de Nantes (LS2N) UMR CNRS 6004
 *
 * @date 2017
 *
 * LUTBasedNSDistanceTransform computes the 2D translated neighborhood-sequence
 * distance transform of a binary image. It reads the input images from its
 * standard input and writes the result to its standard output.
 *
 * This file is part of the DGtal library.
 */

#include <assert.h>
#include <iostream>
#include "Granulometry.h"

/**
 *
 */
Granulometry::Granulometry() : _cols( 0 )
{
}

void
Granulometry::beginOfImage( int cols, int rows )
{
  _cols = cols;
  _maxVal = 0;
  _prevRow =
  (GrayscalePixelType *)malloc( sizeof( GrayscalePixelType ) * cols );
  bzero( _stats, sizeof( _stats ) );
  bzero( _prevRow, cols * sizeof( GrayscalePixelType ) );
}

void
Granulometry::processRow( const GrayscalePixelType * inputRow )
{
  int pix[ 4 ];

  for ( int k = -1; k < _cols; k++ )
  {
    int indices[ 4 ] = {0, 1, 2, 3};
    //
    // | 1 | 2 |
    // | 4 | 8 |
    //
    pix[ 0 ] = k >= 0 ? _prevRow[ k ] : 0;
    pix[ 1 ] = k + 1 < _cols ? _prevRow[ k + 1 ] : 0;
    pix[ 2 ] = k >= 0 && inputRow != NULL ? inputRow[ k ] : 0;
    pix[ 3 ] = k + 1 < _cols && inputRow != NULL ? inputRow[ k + 1 ] : 0;

    _maxVal = pix[ 0 ] > _maxVal ? pix[ 0 ] : _maxVal;

    // Sort indices
    for ( int i = 1; i < 4; i++ )
    {
      for ( int j = 0; j < i; j++ )
      {
        if ( pix[ indices[ j ] ] > pix[ indices[ i ] ] )
        {
          int t = indices[ j ];
          indices[ j ] = indices[ i ];
          indices[ i ] = t;
        }
      }
    }

    int config = 15;
    int threshold = 0;
    for ( int t = 0; t < 4; t++ )
    {
      for ( int i = 0; i < 4; i++ )
      {
        _stats[ i ][ threshold ] += confStats[ config ][ i ];
      }
      threshold = pix[ indices[ t ] ] + 1;
      assert( threshold >= 0 );
      assert( threshold < sizeof( _stats[ 0 ] ) / sizeof( _stats[ 0 ][ 0 ] ) );
      for ( int i = 0; i < 4; i++ )
      {
        _stats[ i ][ threshold ] -= confStats[ config ][ i ];
      }
      config -= 1 << indices[ t ];
    }
    assert( config == 0 );
  }
  // Keep current line
  if ( inputRow != NULL )
    memcpy( _prevRow, inputRow, _cols * sizeof( GrayscalePixelType ) );
  else
    bzero( _prevRow, _cols * sizeof( GrayscalePixelType ) );
}

void
Granulometry::endOfImage()
{
  processRow( NULL );
  for ( int i = 0; i <= _maxVal + 1; i++ )
  {
    _stats[ 3 ][ i ] /= 4;
  }
  /*
  for (int i = 0; i <= _maxVal+1; i++) {
      for (int j = 0; j < 4; j++) {
      printf("%5d\t", _stats[j][i]);
      }
      printf("\n");
  }*/
  for ( int i = 0; i <= _maxVal + 1; i++ )
  {
    for ( int j = 0; j < 4; j++ )
    {
      printf( "%5d\t", _stats[ j ][ i ] );
      _stats[ j ][ i + 1 ] += _stats[ j ][ i ];
    }
    printf( "\n" );
  }
  free( _prevRow );
}

int Granulometry::confStats[ 16 ][ 4 ] = {
{0, 0, 0, 0},  //
{1, 0, 0, 1},  // ▘
{0, 0, 0, 1},  // ▝
{1, 1, 0, 0},  // ▀
{0, 0, 0, 1},  // ▖
{1, 1, 0, 0},  // ▌
{0, 0, 2, -2}, // ▞
{1, 0, 1, -1}, // ▛
{0, 0, 0, 1},  // ▗
{1, 0, 2, -2}, // ▚
{0, 1, 0, 0},  // ▐
{1, 0, 1, -1}, // ▜
{0, 1, 0, 0},  // ▄
{1, 0, 1, -1}, // ▙
{0, 0, 1, -1}, // ▟
{1, 0, 0, 0},  // █
};
/* Old (v1) parameters
int Granulometry::confStats[16][4] =
{
    {0, 0, 0, 0},	//
    {1, 0, 1, 0},	// ▘
    {0, 0, 1, 0},	// ▝
    {1, 1, 0, 0},	// ▀
    {0, 0, 1, 0},	// ▖
    {1, 1, 0, 0},	// ▌
    {0, 0, 0, 2},	// ▞
    {1, 0, 0, 1},	// ▛
    {0, 0, 1, 0},	// ▗
    {1, 0, 0, 2},	// ▚
    {0, 1, 0, 0},	// ▐
    {1, 0, 0, 1},	// ▜
    {0, 1, 0, 0},	// ▄
    {1, 0, 0, 1},	// ▙
    {0, 0, 0, 1},	// ▟
    {1, 0, 0, 0},	// █
};
*/
