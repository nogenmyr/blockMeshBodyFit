/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 Karl-Johan Nogenmyr
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::projectionEdge

Description
    Projects the edge onto a surface using the normals on each point.

SourceFiles
    projectionEdge.C

\*---------------------------------------------------------------------------*/

#ifndef projectionEdge_H
#define projectionEdge_H

#include "curvedEdge.H"
#include "triSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class projectionEdge Declaration
\*---------------------------------------------------------------------------*/

class projectionEdge
:
    public curvedEdge
{
    // Private data

        //- Start and end point. Should have been projected already
        point p1_, p2_;
        
        //- Normals at start and end
        vector n1_, n2_;
        
        //- All the surfaces there are to snap to
        const PtrList<triSurfaceMesh>& snapSurfaces_;
        
        //- How long to search for hit point
        scalar searchLength_;
        
        //- Which surface to snap to
        label surfI_;
        
        
        
    // Private Member Functions


        //- Disallow default bitwise copy construct
        projectionEdge(const projectionEdge&);

        //- Disallow default bitwise assignment
        void operator=(const projectionEdge&);

public:

    //- Runtime type information
//    TypeName("projection");


    // Constructors

        //- Construct from components
        projectionEdge
        (
            const pointField& points,
            const label start, const label end,
            const vector& normalStart,
            const vector& normalEnd,
            const PtrList<triSurfaceMesh>& snapSurfaces,
            const scalar searchLength,
            const label surfI
        );


    //- Destructor
    virtual ~projectionEdge(){}


    // Member Functions

        //- Return the point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        point position(const scalar) const;

        //- Return the length of the curve
        scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
