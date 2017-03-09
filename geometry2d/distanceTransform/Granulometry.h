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
 * @file Granulometry.h
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

#include "ImageFilter.h"
#include "NeighborhoodSequenceDistance.h"

/**
 * An ImageConsumer accepts pixels one row at a time in the standard scan order.
 * A call to beginOfImage() initialises the ImageConsumer,
 * followed by as many calls to processRow() as the number of rows in the image
 * and a final call to endOfImage().
 */

/**
 *
 */
class Granulometry : public ImageConsumer<GrayscalePixelType>
{
public:
  Granulometry( NeighborhoodSequenceDistance * dist );
  void
  beginOfImage( int cols, int rows );
  void
  processRow( const GrayscalePixelType * inputRow );
  void
  endOfImage();

protected:
  int _cols;
  int _maxVal;
  NeighborhoodSequenceDistance * _dist;
  int _stats[ 4 ][ 256 + 1 ];
  GrayscalePixelType * _prevRow;
  static int confStats[ 16 ][ 4 ];
};
