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

\*---------------------------------------------------------------------------*/

#include "projectionEdge.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::projectionEdge::projectionEdge
(
    const pointField& points,
    const label start, const label end,
    const vector& normalStart,
    const vector& normalEnd,
    const PtrList<triSurfaceMesh>& snapSurfaces,
    const scalar searchLength,
    const label surfI
)
:
    curvedEdge(points, start, end),
    p1_(points_[start_]),
    p2_(points_[end_]),
    n1_(normalStart),
    n2_(normalEnd),
    snapSurfaces_(snapSurfaces),
    searchLength_(searchLength),
    surfI_(surfI)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::projectionEdge::position(const scalar lambda) const
{
    if (lambda < 0 || lambda > 1)
    {
        FatalErrorIn("projectionEdge::position(const scalar lambda) const")
            << "Parameter out of range, lambda = " << lambda
            << abort(FatalError);
    }

    if (lambda < SMALL)
    {
        return p1_;
    }
    else if (lambda > 1 - SMALL)
    {
        return p2_;
    }
    else
    {
        point p(lambda*(p2_ - p1_) + p1_);
        vector normal(lambda*(n2_ - n1_) + n1_);
        pointField start(2, p);
        pointField end(2);
        end[0] = p+searchLength_*normal;
        end[1] = p-searchLength_*normal;
        List<pointIndexHit> hits(2);
        snapSurfaces_[surfI_].findLine(start, end, hits);
        scalar hit0Dist = 1e15;
        bool hitFound = false;
        if(hits[0].hit())
        {
            hit0Dist = mag(hits[0].hitPoint()-p);
            hitFound = true;
        }
        scalar hit1Dist = 1e15;
        if(hits[1].hit())
        {
            hit1Dist = mag(hits[1].hitPoint()-p);
            hitFound = true;
        }
        if(hitFound)
        {
            if(hit0Dist < hit1Dist)
            {
                return hits[0].hitPoint();
            }
            else
            {
                return hits[1].hitPoint();
            }
        }
        
        return p;
    }
}


Foam::scalar Foam::projectionEdge::length() const
{
    return 1;
}


// ************************************************************************* //
