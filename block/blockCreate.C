/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | Unsupported Contributions for OpenFOAM
 \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
   \\/     M anipulation  | Copyright (C) 2013 Karl-Johan Nogenmyr
-------------------------------------------------------------------------------
2013-06-01 Karl-Johan Nogenmyr
Implementing:
        label vtxLabelNoBoundary(...)
        vector getPoint(...)
        void getProjection(...)

Complete rewrite of createPoints(...)

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

#include "error.H"
#include "block.H"
#include "Lse.h"
#include "Log.h"
#include "AmgSolver.h"
#include "GaussianElimination.h"

unsigned int loglevel_ = WARNING;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::block::vtxLabel(label i, label j, label k) const
{
    return
    (
        i
      + j * (meshDensity().x() + 1)
      + k * (meshDensity().x() + 1) * (meshDensity().y() + 1)
    );
}


Foam::label Foam::block::vtxLabelNoBoundary(label i, label j, label k) const
{
    return max
    (
        (
            (i-1)
          + (j-1) * (meshDensity().x() - 1)
          + (k-1) * (meshDensity().x() - 1) * (meshDensity().y() - 1)
        ),
        0
    );
}


Foam::vector Foam::block::getPoint(label i, label j, label k, bool flat)  const
{
    const point& p000 = blockPoint(0);
    const point& p100 = blockPoint(1);
    const point& p110 = blockPoint(2);
    const point& p010 = blockPoint(3);

    const point& p001 = blockPoint(4);
    const point& p101 = blockPoint(5);
    const point& p111 = blockPoint(6);
    const point& p011 = blockPoint(7);

    // list of edge point and weighting factors
    const List< List<point> >& p = blockEdgePoints();
    const scalarListList& w = blockEdgeWeights();

    // points on edges
    vector edgex1 = p000 + (p100 - p000)*w[0][i];
    vector edgex2 = p010 + (p110 - p010)*w[1][i];
    vector edgex3 = p011 + (p111 - p011)*w[2][i];
    vector edgex4 = p001 + (p101 - p001)*w[3][i];

    vector edgey1 = p000 + (p010 - p000)*w[4][j];
    vector edgey2 = p100 + (p110 - p100)*w[5][j];
    vector edgey3 = p101 + (p111 - p101)*w[6][j];
    vector edgey4 = p001 + (p011 - p001)*w[7][j];

    vector edgez1 = p000 + (p001 - p000)*w[8][k];
    vector edgez2 = p100 + (p101 - p100)*w[9][k];
    vector edgez3 = p110 + (p111 - p110)*w[10][k];
    vector edgez4 = p010 + (p011 - p010)*w[11][k];
    // calculate the importance factors for all edges

    // x-direction
    scalar impx1 =
    (
        (1.0 - w[0][i])*(1.0 - w[4][j])*(1.0 - w[8][k])
      + w[0][i]*(1.0 - w[5][j])*(1.0 - w[9][k])
    );

    scalar impx2 =
    (
        (1.0 - w[1][i])*w[4][j]*(1.0 - w[11][k])
      + w[1][i]*w[5][j]*(1.0 - w[10][k])
    );

    scalar impx3 =
    (
         (1.0 - w[2][i])*w[7][j]*w[11][k]
       + w[2][i]*w[6][j]*w[10][k]
    );

    scalar impx4 =
    (
        (1.0 - w[3][i])*(1.0 - w[7][j])*w[8][k]
      + w[3][i]*(1.0 - w[6][j])*w[9][k]
    );

    scalar magImpx = impx1 + impx2 + impx3 + impx4;
    impx1 /= magImpx;
    impx2 /= magImpx;
    impx3 /= magImpx;
    impx4 /= magImpx;


    // y-direction
    scalar impy1 =
    (
        (1.0 - w[4][j])*(1.0 - w[0][i])*(1.0 - w[8][k])
      + w[4][j]*(1.0 - w[1][i])*(1.0 - w[11][k])
    );

    scalar impy2 =
    (
        (1.0 - w[5][j])*w[0][i]*(1.0 - w[9][k])
      + w[5][j]*w[1][i]*(1.0 - w[10][k])
    );

    scalar impy3 =
    (
        (1.0 - w[6][j])*w[3][i]*w[9][k]
      + w[6][j]*w[2][i]*w[10][k]
    );

    scalar impy4 =
    (
        (1.0 - w[7][j])*(1.0 - w[3][i])*w[8][k]
      + w[7][j]*(1.0 - w[2][i])*w[11][k]
    );

    scalar magImpy = impy1 + impy2 + impy3 + impy4;
    impy1 /= magImpy;
    impy2 /= magImpy;
    impy3 /= magImpy;
    impy4 /= magImpy;


    // z-direction
    scalar impz1 =
    (
        (1.0 - w[8][k])*(1.0 - w[0][i])*(1.0 - w[4][j])
      + w[8][k]*(1.0 - w[3][i])*(1.0 - w[7][j])
    );

    scalar impz2 =
    (
        (1.0 - w[9][k])*w[0][i]*(1.0 - w[5][j])
      + w[9][k]*w[3][i]*(1.0 - w[6][j])
    );

    scalar impz3 =
    (
        (1.0 - w[10][k])*w[1][i]*w[5][j]
      + w[10][k]*w[2][i]*w[6][j]
    );

    scalar impz4 =
    (
        (1.0 - w[11][k])*(1.0 - w[1][i])*w[4][j]
      + w[11][k]*(1.0 - w[2][i])*w[7][j]
    );

    scalar magImpz = impz1 + impz2 + impz3 + impz4;
    impz1 /= magImpz;
    impz2 /= magImpz;
    impz3 /= magImpz;
    impz4 /= magImpz;


    // calculate the correction vectors
    vector corx1 = impx1*(p[0][i] - edgex1);
    vector corx2 = impx2*(p[1][i] - edgex2);
    vector corx3 = impx3*(p[2][i] - edgex3);
    vector corx4 = impx4*(p[3][i] - edgex4);

    vector cory1 = impy1*(p[4][j] - edgey1);
    vector cory2 = impy2*(p[5][j] - edgey2);
    vector cory3 = impy3*(p[6][j] - edgey3);
    vector cory4 = impy4*(p[7][j] - edgey4);

    vector corz1 = impz1*(p[8][k] - edgez1);
    vector corz2 = impz2*(p[9][k] - edgez2);
    vector corz3 = impz3*(p[10][k] - edgez3);
    vector corz4 = impz4*(p[11][k] - edgez4);


    // multiply by the importance factor

    // x-direction
    edgex1 *= impx1;
    edgex2 *= impx2;
    edgex3 *= impx3;
    edgex4 *= impx4;

    // y-direction
    edgey1 *= impy1;
    edgey2 *= impy2;
    edgey3 *= impy3;
    edgey4 *= impy4;

    // z-direction
    edgez1 *= impz1;
    edgez2 *= impz2;
    edgez3 *= impz3;
    edgez4 *= impz4;


    // add the contributions
    vector result =
    (
        edgex1 + edgex2 + edgex3 + edgex4
      + edgey1 + edgey2 + edgey3 + edgey4
      + edgez1 + edgez2 + edgez3 + edgez4
    ) / 3.0;

    if (flat)
    {
        return result;
    }

    result +=
    (
        corx1 + corx2 + corx3 + corx4
      + cory1 + cory2 + cory3 + cory4
      + corz1 + corz2 + corz3 + corz4
    );
    return result;
}

void Foam::block::getProjection
(
    point& p,
    vector& normal,
    const triSurfaceMesh& snapSurface,
    const scalar searchLength
) const
{
    pointField start(2, p);
    pointField end(2);
    end[0] = p + searchLength*normal;
    end[1] = p - searchLength*normal;
    List<pointIndexHit> hits(2);
    snapSurface.findLine(start, end, hits);
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
            p = hits[0].hitPoint();
        }
        else
        {
            p = hits[1].hitPoint();
        }
    }
}

void Foam::block::createPoints
(
    const PtrList<triSurfaceMesh>& snapSurfaces, //the stls
    const labelList& snapList,  //which face should snap to which surface (or not -> =-1)
    const faceList& allBlocksFaces,  // all blocks' faces
    const vectorField& pointNormals, //blocks' corner's projection normals
    const scalar searchLength
) const
{
    // set local variables for mesh specification
    const label ni = meshDensity().x();
    const label nj = meshDensity().y();
    const label nk = meshDensity().z();

    const faceList thisBlocksFaces(blockShape().faces());
    labelList thisBlocksSnaps(thisBlocksFaces.size(), -1);

    forAll(snapList, sfaceI)
    {
        label surfI = snapList[sfaceI];
        if(surfI>=0)
        {
            forAll(thisBlocksFaces, bfaceI)
            {
                bool match = true;
                forAll(allBlocksFaces[sfaceI], point1)
                {
                    label face1p = allBlocksFaces[sfaceI][point1];
                    bool pointFound = false;
                    forAll(thisBlocksFaces[bfaceI], point2)
                    {
                        label face2p = thisBlocksFaces[bfaceI][point2];
                        if(face1p == face2p)
                        {
                            pointFound = true;
                        }
                    }
                    match *= pointFound;
                }
                if (match)
                {
                    thisBlocksSnaps[bfaceI] = surfI;
                }
            }
        }
    }
    
    // list of edge weighting factors
    const scalarListList& w = blockEdgeWeights();

    // Generate surface mesh first:
    //    Mesh of flat faces
    List< List<vector > > i0(nj+1);
    List< List<vector > > i1(nj+1);

    List< List<vector > > j0(ni+1);
    List< List<vector > > j1(ni+1);

    List< List<vector > > k1(ni+1);
    List< List<vector > > k0(ni+1);

    //    ... and mesh of projected faces
    List< List<vector > > pi0(nj+1);
    List< List<vector > > pi1(nj+1);

    List< List<vector > > pj0(ni+1);
    List< List<vector > > pj1(ni+1);

    List< List<vector > > pk1(ni+1);
    List< List<vector > > pk0(ni+1);

    for (label j = 0; j <= nj; j++)
    {
        i0[j].resize(nk+1);
        i1[j].resize(nk+1);
        pi0[j].resize(nk+1);
        pi1[j].resize(nk+1);
        for (label k = 0; k <= nk; k++)
        {
// First plane i0
            point p = getPoint(0,j,k, false);// intp surface from edges only
            
            if(thisBlocksSnaps[0]>=0)
            {
//                p = getPoint(0,j,k, true);
                vector normal
                (
                    pointNormals[blockShape()[0]]*(1-w[4][j])*(1-w[8][k])
                  + pointNormals[blockShape()[3]]*w[4][j]*(1-w[11][k])
                  + pointNormals[blockShape()[7]]*w[11][k]*w[7][j]
                  + pointNormals[blockShape()[4]]*(1-w[7][j])*w[8][k]
                );  
                getProjection
                (
                    p, 
                    normal, 
                    snapSurfaces[thisBlocksSnaps[0]], 
                    searchLength
                );
            }

            pi0[j][k] = p; //final surface point
            
            p = getPoint(ni,j,k, false);// intp surface from edges only

            if(thisBlocksSnaps[1]>=0)
            {
//                p = getPoint(ni,j,k, true);
                vector normal
                (
                    pointNormals[blockShape()[1]]*(1-w[5][j])*(1-w[9][k])
                  + pointNormals[blockShape()[2]]*w[5][j]*(1-w[10][k])
                  + pointNormals[blockShape()[6]]*w[10][k]*w[6][j]
                  + pointNormals[blockShape()[5]]*(1-w[6][j])*w[9][k]
                );  
                getProjection
                (
                    p, 
                    normal, 
                    snapSurfaces[thisBlocksSnaps[1]],
                    searchLength
                );
            }

            pi1[j][k] = p;

            i0[j][k] = getPoint(0,j,k, true);//flat
            i1[j][k] = getPoint(ni,j,k, true);//flat
        }
    }

    for (label i = 0; i <= ni; i++)
    {
        j0[i].resize(nk+1);
        j1[i].resize(nk+1);
        pj0[i].resize(nk+1);
        pj1[i].resize(nk+1);
        for (label k = 0; k <= nk; k++)
        {
            point p = getPoint(i,0,k, false);

            if(thisBlocksSnaps[2]>=0)
            {
//                p = getPoint(i,0,k, true);
                vector normal
                (
                    pointNormals[blockShape()[0]]*(1-w[0][i])*(1-w[8][k])
                  + pointNormals[blockShape()[1]]*w[0][i]*(1-w[9][k])
                  + pointNormals[blockShape()[5]]*w[9][k]*w[3][i]
                  + pointNormals[blockShape()[4]]*(1-w[3][i])*w[8][k]
                );  
                getProjection
                (
                    p, 
                    normal, 
                    snapSurfaces[thisBlocksSnaps[2]],
                    searchLength
                );
            }

            pj0[i][k] = p;
            
            p = getPoint(i,nj,k, false);

            if(thisBlocksSnaps[3]>=0)
            {
//                p = getPoint(i,nj,k, true);
                vector normal
                (
                    pointNormals[blockShape()[3]]*(1-w[1][i])*(1-w[11][k])
                  + pointNormals[blockShape()[2]]*w[1][i]*(1-w[10][k])
                  + pointNormals[blockShape()[6]]*w[2][i]*w[10][k]
                  + pointNormals[blockShape()[7]]*(1-w[2][i])*w[11][k]
                );  
                getProjection
                (
                    p, 
                    normal, 
                    snapSurfaces[thisBlocksSnaps[3]],
                    searchLength
                );
            }

            pj1[i][k] = p;

            j0[i][k] = getPoint(i,0,k, true);
            j1[i][k] = getPoint(i,nj,k, true);
        }
    }


    for (label i = 0; i <= ni; i++)
    {
        k0[i].resize(nj+1);
        k1[i].resize(nj+1);
        pk0[i].resize(nj+1);
        pk1[i].resize(nj+1);
        for (label j = 0; j <= nj; j++)
        {
            point p = getPoint(i,j,0, false);

            if(thisBlocksSnaps[4]>=0)
            {
//                p = getPoint(i,j,0, true);
                vector normal
                (
                    pointNormals[blockShape()[0]]*(1-w[0][i])*(1-w[4][j])
                  + pointNormals[blockShape()[1]]*w[0][i]*(1-w[5][j])
                  + pointNormals[blockShape()[2]]*w[5][j]*w[1][i]
                  + pointNormals[blockShape()[3]]*(1-w[1][i])*w[4][j]
                );  
                getProjection
                (
                    p, 
                    normal, 
                    snapSurfaces[thisBlocksSnaps[4]],
                    searchLength
                );
            }
            pk0[i][j] = p;

            p = getPoint(i,j,nk, false);

            if(thisBlocksSnaps[5]>=0)
            {
//                p = getPoint(i,j,nk, true);
                vector normal
                (
                    pointNormals[blockShape()[4]]*(1-w[3][i])*(1-w[7][j])
                  + pointNormals[blockShape()[5]]*w[3][i]*(1-w[6][j])
                  + pointNormals[blockShape()[6]]*w[6][j]*w[2][i]
                  + pointNormals[blockShape()[7]]*(1-w[2][i])*w[7][j]
                );  
                getProjection
                (
                    p, 
                    normal, 
                    snapSurfaces[thisBlocksSnaps[5]],
                    searchLength
                );
            }
            pk1[i][j] = p;

            k0[i][j] = getPoint(i,j,0, true);
            k1[i][j] = getPoint(i,j,nk, true);
        }
    }

    //
    // generate vertices
    //
    vertices_.clear();
    vertices_.setSize(nPoints());
    
    // Set up the generic multigrid solver for solving Poisson eqs
	Lse linSystemX;
	Lse linSystemY;
	Lse linSystemZ;
	
	label N = (ni-1)*(nj-1)*(nk-1);

    std::vector<std::size_t> capacities(N, 19);
	linSystemX.b_.resize(N, false);
    linSystemX.x_.resize(N, false);

	linSystemY.b_.resize(N, false);
    linSystemY.x_.resize(N, false);

	linSystemZ.b_.resize(N, false);
    linSystemZ.x_.resize(N, false);

    Info << "Solving system of size " << N <<endl;

    // First place the point algebraically
    for (label k = 0; k <= nk; k++)
    {
        for (label j = 0; j <= nj; j++)
        {
            for (label i = 0; i <= ni; i++)
            {
                const label vertexNo = vtxLabel(i, j, k);
                if (k==0 || k==nk)
                {
                    vertices_[vertexNo] = (1-w[8][k])*pk0[i][j] + w[8][k]*pk1[i][j];
                }
                else if (j==0 || j==nj)
                {
                    vertices_[vertexNo] = (1-w[4][j])*pj0[i][k] + w[4][j]*pj1[i][k];
                }
                else if (i==0 || i==ni)
                {
                    vertices_[vertexNo] = (1-w[0][i])*pi0[j][k] + w[0][i]*pi1[j][k];
                }
                else
                {
                    scalar g = (1-w[8][k])*(1-w[4][j]) + w[4][j]*(1-w[9][k]) + w[7][j]*w[9][k] + w[8][k]*(1-w[7][j]);
                    scalar f =  (1-w[11][k])*(1-w[5][j]) + w[5][j]*(1-w[10][k]) + w[6][j]*w[10][k] + w[11][k]*(1-w[6][j]);
                    
                    scalar y = (1-w[8][k])*(1-w[4][j])*w[0][i] + w[4][j]*(1-w[9][k])*w[1][i]
                             + w[7][j]*w[9][k]*w[2][i] + w[8][k]*(1-w[7][j])*w[3][i];
                             
                    scalar x = (1-w[11][k])*(1-w[5][j])*w[0][i] + w[5][j]*(1-w[10][k])*w[1][i] +
                               w[6][j]*w[10][k]*w[2][i] + w[11][k]*(1-w[6][j])*w[3][i];
                    
                    scalar wi = 0.;
                    if (abs(f-g)>SMALL)
                    {
                        wi = 0.5/(f-g) * (-g+x-y + sqrt(pow(-g+x-y,2) + 4.*(f-g)*y));
                    }
                    else
                    {
                        wi = y /(g-x+y);
                    }

                    g = (1-w[0][i])*(1-w[8][k]) + w[0][i]*(1-w[9][k]) + w[3][i]*w[9][k] + w[8][k]*(1-w[3][i]);
                    
                    f = (1-w[11][k])*(1-w[1][i]) + w[1][i]*(1-w[10][k]) + w[2][i]*w[10][k] + w[11][k]*(1-w[2][i]);
                    
                    y = (1-w[0][i])*(1-w[8][k])*w[4][j] + w[0][i]*(1-w[9][k])*w[5][j] + w[3][i]*w[9][k]*w[6][j] + w[8][k]*(1-w[3][i])*w[7][j];
                    
                    x = (1-w[11][k])*(1-w[1][i])*w[4][j] + w[1][i]*(1-w[10][k])*w[5][j] + w[2][i]*w[10][k]*w[6][j] + w[11][k]*(1-w[2][i])*w[7][j];
                    
                    scalar wj = 0.;
                    if (abs(f-g)>SMALL)
                    {
                        wj = 0.5/(f-g) * (-g+x-y + sqrt(pow(-g+x-y,2) + 4.*(f-g)*y));
                    }
                    else
                    {
                        wj = y /(g-x+y);
                    }

                    g = (1-w[4][j])*(1-w[0][i]) + w[0][i]*(1-w[5][j]) + w[1][i]*w[5][j] + w[4][j]*(1-w[1][i]);
                    
                    f = (1-w[7][j])*(1-w[3][i]) + w[3][i]*(1-w[6][j]) + w[2][i]*w[6][j] + w[7][j]*(1-w[2][i]);
                    
                    y = (1-w[4][j])*(1-w[0][i])*w[8][k] + w[0][i]*(1-w[5][j])*w[9][k] + w[1][i]*w[5][j]*w[10][k] + w[4][j]*(1-w[1][i])*w[11][k];
                    
                    x = (1-w[7][j])*(1-w[3][i])*w[8][k] + w[3][i]*(1-w[6][j])*w[9][k] + w[2][i]*w[6][j]*w[10][k] + w[7][j]*(1-w[2][i])*w[11][k];
                    
                    scalar wk = 0.;
                    if (abs(f-g)>SMALL)
                    {
                        wk = 0.5/(f-g) * (-g+x-y + sqrt(pow(-g+x-y,2) + 4.*(f-g)*y));
                    }
                    else
                    {
                        wk = y /(g-x+y);
                    }

                    vector iline = (1-wi)*i0[j][k] + wi*i1[j][k];
                    vector jline = (1-wj)*j0[i][k] + wj*j1[i][k];
                    vector kline = (1-wk)*k0[i][j] + wk*k1[i][j];

                    const label vertexNo = vtxLabel(i, j, k);
                    
                    vertices_[vertexNo] = 
                    (
                         wj*(1-wj)*wk*(1-wk)*iline
                       + wi*(1-wi)*wk*(1-wk)*jline
                       + wi*(1-wi)*wj*(1-wj)*kline
                    )/
                    (
                         wj*(1-wj)*wk*(1-wk)
                       + wi*(1-wi)*wk*(1-wk) 
                       + wi*(1-wi)*wj*(1-wj)
                    );

                    vector piline = (1-wi)*(pi0[j][k]-i0[j][k]) + wi*(pi1[j][k]-i1[j][k]);
                    vector pjline = (1-wj)*(pj0[i][k]-j0[i][k]) + wj*(pj1[i][k]-j1[i][k]);
                    vector pkline = (1-wk)*(pk0[i][j]-k0[i][j]) + wk*(pk1[i][j]-k1[i][j]);


                    vertices_[vertexNo] += 
                    (
                         wj*(1-wj)*wk*(1-wk)*piline
                       + wi*(1-wi)*wk*(1-wk)*pjline
                       + wi*(1-wi)*wj*(1-wj)*pkline
                    )/
                    (
                         wj*(1-wj)*wk*(1-wk)
                       + wi*(1-wi)*wk*(1-wk) 
                       + wi*(1-wi)*wj*(1-wj)
                    );
                    
                    const label nbijk = vtxLabelNoBoundary(i, j, k);

                    linSystemX.x_[nbijk] = vertices_[vertexNo][0];
                    linSystemY.x_[nbijk] = vertices_[vertexNo][1];
                    linSystemZ.x_[nbijk] = vertices_[vertexNo][2];
                }
            }
        }
    }

// Now we have placed the points algebraically. This is fine for boundaries of each block, but
// generally not good enough for block's interior points. For this we need to solve a PDE.
// The algebraically determined point positions are used as initial and boundary conditions.

// Following Thompson, Joe F., Handbook of Grid Generation (1999) Chapter 4.4

    if (ni == 1 || nj == 1 || nk == 1) return; // no interior points exists!
    
    scalar dxi = 1./(ni-1.);
    scalar deta = 1./(nj-1.);
    scalar dzeta = 1./(nk-1.);

    vectorField grading(nPoints()); // Parameter space P
    for (label k = 0; k <= nk; k++)
    {
        for (label j = 0; j <= nj; j++)
        {
            for (label i = 0; i <= ni; i++)
            {
                const label ijk = vtxLabel(i, j, k);
                const point p000(0,0,0);
                const point p100(1,0,0);
                const point p110(1,1,0);
                const point p010(0,1,0);

                const point p001(0,0,1);
                const point p101(1,0,1);
                const point p111(1,1,1);
                const point p011(0,1,1);

                // list of edge weighting factors
                const scalarListList& w = blockEdgeWeights();

                // points on edges
                vector edgex1 = p000 + (p100 - p000)*w[0][i];
                vector edgex2 = p010 + (p110 - p010)*w[1][i];
                vector edgex3 = p011 + (p111 - p011)*w[2][i];
                vector edgex4 = p001 + (p101 - p001)*w[3][i];

                vector edgey1 = p000 + (p010 - p000)*w[4][j];
                vector edgey2 = p100 + (p110 - p100)*w[5][j];
                vector edgey3 = p101 + (p111 - p101)*w[6][j];
                vector edgey4 = p001 + (p011 - p001)*w[7][j];

                vector edgez1 = p000 + (p001 - p000)*w[8][k];
                vector edgez2 = p100 + (p101 - p100)*w[9][k];
                vector edgez3 = p110 + (p111 - p110)*w[10][k];
                vector edgez4 = p010 + (p011 - p010)*w[11][k];
                // calculate the importance factors for all edges

                // x-direction
                scalar impx1 =
                (
                    (1.0 - w[0][i])*(1.0 - w[4][j])*(1.0 - w[8][k])
                  + w[0][i]*(1.0 - w[5][j])*(1.0 - w[9][k])
                );

                scalar impx2 =
                (
                    (1.0 - w[1][i])*w[4][j]*(1.0 - w[11][k])
                  + w[1][i]*w[5][j]*(1.0 - w[10][k])
                );

                scalar impx3 =
                (
                     (1.0 - w[2][i])*w[7][j]*w[11][k]
                   + w[2][i]*w[6][j]*w[10][k]
                );

                scalar impx4 =
                (
                    (1.0 - w[3][i])*(1.0 - w[7][j])*w[8][k]
                  + w[3][i]*(1.0 - w[6][j])*w[9][k]
                );

                scalar magImpx = impx1 + impx2 + impx3 + impx4;
                impx1 /= magImpx;
                impx2 /= magImpx;
                impx3 /= magImpx;
                impx4 /= magImpx;


                // y-direction
                scalar impy1 =
                (
                    (1.0 - w[4][j])*(1.0 - w[0][i])*(1.0 - w[8][k])
                  + w[4][j]*(1.0 - w[1][i])*(1.0 - w[11][k])
                );

                scalar impy2 =
                (
                    (1.0 - w[5][j])*w[0][i]*(1.0 - w[9][k])
                  + w[5][j]*w[1][i]*(1.0 - w[10][k])
                );

                scalar impy3 =
                (
                    (1.0 - w[6][j])*w[3][i]*w[9][k]
                  + w[6][j]*w[2][i]*w[10][k]
                );

                scalar impy4 =
                (
                    (1.0 - w[7][j])*(1.0 - w[3][i])*w[8][k]
                  + w[7][j]*(1.0 - w[2][i])*w[11][k]
                );

                scalar magImpy = impy1 + impy2 + impy3 + impy4;
                impy1 /= magImpy;
                impy2 /= magImpy;
                impy3 /= magImpy;
                impy4 /= magImpy;


                // z-direction
                scalar impz1 =
                (
                    (1.0 - w[8][k])*(1.0 - w[0][i])*(1.0 - w[4][j])
                  + w[8][k]*(1.0 - w[3][i])*(1.0 - w[7][j])
                );

                scalar impz2 =
                (
                    (1.0 - w[9][k])*w[0][i]*(1.0 - w[5][j])
                  + w[9][k]*w[3][i]*(1.0 - w[6][j])
                );

                scalar impz3 =
                (
                    (1.0 - w[10][k])*w[1][i]*w[5][j]
                  + w[10][k]*w[2][i]*w[6][j]
                );

                scalar impz4 =
                (
                    (1.0 - w[11][k])*(1.0 - w[1][i])*w[4][j]
                  + w[11][k]*(1.0 - w[2][i])*w[7][j]
                );

                scalar magImpz = impz1 + impz2 + impz3 + impz4;
                impz1 /= magImpz;
                impz2 /= magImpz;
                impz3 /= magImpz;
                impz4 /= magImpz;

                // multiply by the importance factor

                // x-direction
                edgex1 *= impx1;
                edgex2 *= impx2;
                edgex3 *= impx3;
                edgex4 *= impx4;

                // y-direction
                edgey1 *= impy1;
                edgey2 *= impy2;
                edgey3 *= impy3;
                edgey4 *= impy4;

                // z-direction
                edgez1 *= impz1;
                edgez2 *= impz2;
                edgez3 *= impz3;
                edgez4 *= impz4;


                // add the contributions
                grading[ijk] =
                (
                    edgex1 + edgex2 + edgex3 + edgex4
                  + edgey1 + edgey2 + edgey3 + edgey4
                  + edgez1 + edgez2 + edgez3 + edgez4
                ) / 3.0;
            }
        }
    }
    
    // Grid control functions
    vectorField P11(nPoints());
    vectorField P12(nPoints());
    vectorField P13(nPoints());
    vectorField P22(nPoints());
    vectorField P23(nPoints());
    vectorField P33(nPoints());

    for (label k = 1; k < nk; k++)
    {
        for (label j = 1; j < nj; j++)
        {
            for (label i = 1; i < ni; i++)
            {
                const label ijk = vtxLabel(i, j, k);

                const label ip1 = vtxLabel(i+1, j, k); //forward
                const label jp1 = vtxLabel(i, j+1, k);
                const label kp1 = vtxLabel(i, j, k+1);

                const label im1 = vtxLabel(i-1, j, k); //backward
                const label jm1 = vtxLabel(i, j-1, k);
                const label km1 = vtxLabel(i, j, k-1);

                const label im1jm1 = vtxLabel(i-1, j-1, k);
                const label im1jp1 = vtxLabel(i-1, j+1, k);
                const label ip1jm1 = vtxLabel(i+1, j-1, k);
                const label ip1jp1 = vtxLabel(i+1, j+1, k);

                const label im1km1 = vtxLabel(i-1, j, k-1);
                const label im1kp1 = vtxLabel(i-1, j, k+1);
                const label ip1km1 = vtxLabel(i+1, j, k-1);
                const label ip1kp1 = vtxLabel(i+1, j, k+1);

                const label jm1km1 = vtxLabel(i, j-1, k-1);
                const label jm1kp1 = vtxLabel(i, j-1, k+1);
                const label jp1km1 = vtxLabel(i, j+1, k-1);
                const label jp1kp1 = vtxLabel(i, j+1, k+1);
                
                vector dstudxi =   (grading[ip1] - grading[im1] )/(2.*dxi);
                vector dstudeta =  (grading[jp1] - grading[jm1] )/(2.*deta);
                vector dstudzeta = (grading[kp1] - grading[km1] )/(2.*dzeta);

                vector d2studxi2 =   (grading[ip1] - 2*grading[ijk] + grading[im1] )/(dxi*dxi);
                vector d2studeta2 =  (grading[jp1] - 2*grading[ijk] + grading[jm1] )/(deta*deta);
                vector d2studzeta2 = (grading[kp1] - 2*grading[ijk] + grading[km1] )/(dzeta*dzeta);

                vector d2studxideta = (grading[im1jm1]  - grading[im1jp1] - grading[ip1jm1] + grading[ip1jp1])/(4*deta*dxi);
                vector d2studxidzeta = (grading[im1km1]  - grading[im1kp1] - grading[ip1km1] + grading[ip1kp1])/(4*dzeta*dxi);
                vector d2studetadzeta = (grading[jm1km1]  - grading[jm1kp1] - grading[jp1km1] + grading[jp1kp1])/(4*dzeta*deta);
                
                tensor T
                (
                    dstudxi[0], dstudeta[0], dstudzeta[0],
                    dstudxi[1], dstudeta[1], dstudzeta[1], 
                    dstudxi[2], dstudeta[2], dstudzeta[2]
                );

                tensor Tinv = inv(T);

                P11[ijk] = -Tinv & d2studxi2;
                P22[ijk] = -Tinv & d2studeta2;
                P33[ijk] = -Tinv & d2studzeta2;
                
                P12[ijk] = -Tinv & d2studxideta;
                P13[ijk] = -Tinv & d2studxidzeta;
                P23[ijk] = -Tinv & d2studetadzeta;
            }
        }
    }

    Info <<"Solving Poisson";
    for (int gridIters=0; gridIters<2; gridIters++)
    {
        Info << ".";
	    int nu1 = 1, nu2 = 1;
	    scalar thetaPos = 0.1, thetaNeg = 0.4;
        // setup AMG solver
        AmgSolver<GaussianElimination> amg;

        // Full Multigrid
        amg.setStartupExecution(false);

        // V(1,1)-cycle
        amg.setSmoothingParameters(nu1, nu2);

        // set strong dependency thresholds for positive and negative entries
        amg.setCoarseningParameters(thetaPos, thetaNeg);

        // choose coarsening strategy from coarsening_amg1r5, coarsening_direct and coarsening_standard
        amg.setCoarseningStrategy(AmgSolver<GaussianElimination>::coarsening_direct);
        
        linSystemX.A_ = MatrixCRS(N, N, capacities);
        for (label k = 1; k < nk; k++)
        {
            for (label j = 1; j < nj; j++)
            {
                for (label i = 1; i < ni; i++)
                {
                    const label ijk = vtxLabel(i, j, k);

                    const label ip1 = vtxLabel(i+1, j, k); //forward
                    const label jp1 = vtxLabel(i, j+1, k);
                    const label kp1 = vtxLabel(i, j, k+1);

                    const label im1 = vtxLabel(i-1, j, k); //backward
                    const label jm1 = vtxLabel(i, j-1, k);
                    const label km1 = vtxLabel(i, j, k-1);

                    const label im1jm1 = vtxLabel(i-1, j-1, k);
                    const label im1jp1 = vtxLabel(i-1, j+1, k);
                    const label ip1jm1 = vtxLabel(i+1, j-1, k);
                    const label ip1jp1 = vtxLabel(i+1, j+1, k);

                    const label im1km1 = vtxLabel(i-1, j, k-1);
                    const label im1kp1 = vtxLabel(i-1, j, k+1);
                    const label ip1km1 = vtxLabel(i+1, j, k-1);
                    const label ip1kp1 = vtxLabel(i+1, j, k+1);

                    const label jm1km1 = vtxLabel(i, j-1, k-1);
                    const label jm1kp1 = vtxLabel(i, j-1, k+1);
                    const label jp1km1 = vtxLabel(i, j+1, k-1);
                    const label jp1kp1 = vtxLabel(i, j+1, k+1);
                    
                    label nbijk = vtxLabelNoBoundary(i, j, k);//point
                    
                    label nbip1 = vtxLabelNoBoundary(i+1, j, k); //forward
                    label nbjp1 = vtxLabelNoBoundary(i, j+1, k);
                    label nbkp1 = vtxLabelNoBoundary(i, j, k+1);

                    label nbim1 = vtxLabelNoBoundary(i-1, j, k); //backward
                    label nbjm1 = vtxLabelNoBoundary(i, j-1, k);
                    label nbkm1 = vtxLabelNoBoundary(i, j, k-1);

                    label nbim1jm1 = vtxLabelNoBoundary(i-1, j-1, k);
                    label nbim1jp1 = vtxLabelNoBoundary(i-1, j+1, k);
                    label nbip1jm1 = vtxLabelNoBoundary(i+1, j-1, k);
                    label nbip1jp1 = vtxLabelNoBoundary(i+1, j+1, k);

                    label nbim1km1 = vtxLabelNoBoundary(i-1, j, k-1);
                    label nbim1kp1 = vtxLabelNoBoundary(i-1, j, k+1);
                    label nbip1km1 = vtxLabelNoBoundary(i+1, j, k-1);
                    label nbip1kp1 = vtxLabelNoBoundary(i+1, j, k+1);

                    label nbjm1km1 = vtxLabelNoBoundary(i, j-1, k-1);
                    label nbjm1kp1 = vtxLabelNoBoundary(i, j-1, k+1);
                    label nbjp1km1 = vtxLabelNoBoundary(i, j+1, k-1);
                    label nbjp1kp1 = vtxLabelNoBoundary(i, j+1, k+1);

                    vector dxdxi = vector::zero;
                    vector dxdeta = vector::zero;
                    vector dxdzeta = vector::zero;
                    // First derivatives
                    if (i>1 && i<ni-1)
                    {
                        dxdxi = vector
                        (
                            (linSystemX.x_[nbip1] - linSystemX.x_[nbim1])/(2.*dxi), 
                            (linSystemY.x_[nbip1] - linSystemY.x_[nbim1])/(2.*dxi), 
                            (linSystemZ.x_[nbip1] - linSystemZ.x_[nbim1])/(2.*dxi)
                        );
                    }
                    else if (i==1)
                    {
                        dxdxi = vector
                        (
                            (linSystemX.x_[nbip1] - vertices_[im1][0])/(2.*dxi), 
                            (linSystemY.x_[nbip1] - vertices_[im1][1])/(2.*dxi), 
                            (linSystemZ.x_[nbip1] - vertices_[im1][2])/(2.*dxi)
                        );
                    }
                    else
                    {
                        dxdxi = vector
                        (
                            (vertices_[ip1][0] - linSystemX.x_[nbim1])/(2.*dxi), 
                            (vertices_[ip1][1] - linSystemY.x_[nbim1])/(2.*dxi), 
                            (vertices_[ip1][2] - linSystemZ.x_[nbim1])/(2.*dxi)
                        );
                    }

                    if (j>1 && j<nj-1)
                    {
                        dxdeta = vector
                        (
                            (linSystemX.x_[nbjp1] - linSystemX.x_[nbjm1])/(2.*deta),
                            (linSystemY.x_[nbjp1] - linSystemY.x_[nbjm1])/(2.*deta),
                            (linSystemZ.x_[nbjp1] - linSystemZ.x_[nbjm1])/(2.*deta)
                        );
                    }
                    else if (j==1)
                    {
                        dxdeta = vector
                        (
                            (linSystemX.x_[nbjp1] - vertices_[jm1][0])/(2.*deta),
                            (linSystemY.x_[nbjp1] - vertices_[jm1][1])/(2.*deta),
                            (linSystemZ.x_[nbjp1] - vertices_[jm1][2])/(2.*deta)
                        );
                    }
                    else
                    {
                        dxdeta = vector
                        (
                            (vertices_[jp1][0] - linSystemX.x_[nbjm1])/(2.*deta),
                            (vertices_[jp1][1] - linSystemY.x_[nbjm1])/(2.*deta),
                            (vertices_[jp1][2] - linSystemZ.x_[nbjm1])/(2.*deta)
                        );
                    }

                    if (k>1 && k<nk-1)
                    {
                        dxdzeta = vector
                        (
                            (linSystemX.x_[nbkp1] - linSystemX.x_[nbkm1])/(2.*dzeta),
                            (linSystemY.x_[nbkp1] - linSystemY.x_[nbkm1])/(2.*dzeta),
                            (linSystemZ.x_[nbkp1] - linSystemZ.x_[nbkm1])/(2.*dzeta)
                        );
                    }
                    else if (k==1)
                    {
                        dxdzeta = vector
                        (
                            (linSystemX.x_[nbkp1] - vertices_[km1][0])/(2.*dzeta),
                            (linSystemY.x_[nbkp1] - vertices_[km1][1])/(2.*dzeta),
                            (linSystemZ.x_[nbkp1] - vertices_[km1][2])/(2.*dzeta)
                        );
                    }
                    else
                    {
                        dxdzeta = vector
                        (
                            (vertices_[kp1][0] - linSystemX.x_[nbkm1])/(2.*dzeta),
                            (vertices_[kp1][1] - linSystemY.x_[nbkm1])/(2.*dzeta),
                            (vertices_[kp1][2] - linSystemZ.x_[nbkm1])/(2.*dzeta)
                        );
                    }

                    // covariant tensor coeffs
                    scalar a11cov = dxdxi & dxdxi;
                    scalar a12cov = dxdxi & dxdeta;
                    scalar a13cov = dxdxi & dxdzeta;
                    scalar a22cov = dxdeta & dxdeta;
                    scalar a23cov = dxdeta & dxdzeta;
                    scalar a33cov = dxdzeta & dxdzeta;

                    // contravariant tensor coeffs
                    scalar a11 = a22cov*a33cov - a23cov*a23cov;
                    scalar a22 = a11cov*a33cov - a13cov*a13cov;
                    scalar a33 = a11cov*a22cov - a12cov*a12cov;

                    scalar a12 = a13cov*a23cov - a12cov*a33cov;
                    scalar a13 = a12cov*a23cov - a13cov*a22cov;
                    scalar a23 = a13cov*a12cov - a11cov*a23cov;
                    
                    // Grid control
                    vector gridCon = 
                    (
                        a11*P11[ijk]
                      + 2*a12*P12[ijk] 
                      + 2*a13*P13[ijk] 
                      + a22*P22[ijk] 
                      + 2*a23*P23[ijk] 
                      +a33*P33[ijk]
                    );

                    // Matrix coefficients A...
                    scalar Aijk = 
                    (
                      - 2*a11/(dxi*dxi) 
                      - 2*a22/(deta*deta) 
                      - 2*a33/(dzeta*dzeta)
                    );

                    scalar Aip1 = a11/(dxi*dxi) + gridCon[0]*0.5/dxi;
                    scalar Aim1 = a11/(dxi*dxi) - gridCon[0]*0.5/dxi;
                    
                    scalar Ajp1 = a22/(deta*deta) + gridCon[1]*0.5/deta;
                    scalar Ajm1 = a22/(deta*deta) - gridCon[1]*0.5/deta;
                    
                    scalar Akp1 = a33/(dzeta*dzeta) + gridCon[2]*0.5/dzeta;
                    scalar Akm1 = a33/(dzeta*dzeta) - gridCon[2]*0.5/dzeta;
                    
                    scalar Aim1jm1 = +a12/(2*deta*dxi);
                    scalar Aim1jp1 = -a12/(2*deta*dxi);
                    scalar Aip1jm1 = -a12/(2*deta*dxi);
                    scalar Aip1jp1 = +a12/(2*deta*dxi);

                    scalar Aim1km1 = +a13/(2*dzeta*dxi);
                    scalar Aim1kp1 = -a13/(2*dzeta*dxi);
                    scalar Aip1km1 = -a13/(2*dzeta*dxi);
                    scalar Aip1kp1 = +a13/(2*dzeta*dxi);

                    scalar Ajm1km1 = +a23/(2*dzeta*deta);
                    scalar Ajm1kp1 = -a23/(2*dzeta*deta);
                    scalar Ajp1km1 = -a23/(2*dzeta*deta);
                    scalar Ajp1kp1 = +a23/(2*dzeta*deta);

                    //Reset RHS
                    linSystemX.b_[nbijk] = 0.;
                    linSystemY.b_[nbijk] = 0.;
                    linSystemZ.b_[nbijk] = 0.;
                    

                    if (i==1)//Boundary points
                    {
                        linSystemX.b_[nbijk] -= vertices_[im1][0]*Aim1;
                        linSystemY.b_[nbijk] -= vertices_[im1][1]*Aim1;
                        linSystemZ.b_[nbijk] -= vertices_[im1][2]*Aim1;
                        nbim1 = -1;
                        linSystemX.b_[nbijk] -= vertices_[im1jm1][0]*Aim1jm1;
                        linSystemY.b_[nbijk] -= vertices_[im1jm1][1]*Aim1jm1;
                        linSystemZ.b_[nbijk] -= vertices_[im1jm1][2]*Aim1jm1;
                        nbim1jm1 = -1;
                        linSystemX.b_[nbijk] -= vertices_[im1jp1][0]*Aim1jp1;
                        linSystemY.b_[nbijk] -= vertices_[im1jp1][1]*Aim1jp1;
                        linSystemZ.b_[nbijk] -= vertices_[im1jp1][2]*Aim1jp1;
                        nbim1jp1 = -1;
                        linSystemX.b_[nbijk] -= vertices_[im1km1][0]*Aim1km1;
                        linSystemY.b_[nbijk] -= vertices_[im1km1][1]*Aim1km1;
                        linSystemZ.b_[nbijk] -= vertices_[im1km1][2]*Aim1km1;
                        nbim1km1 = -1;
                        linSystemX.b_[nbijk] -= vertices_[im1kp1][0]*Aim1kp1;
                        linSystemY.b_[nbijk] -= vertices_[im1kp1][1]*Aim1kp1;
                        linSystemZ.b_[nbijk] -= vertices_[im1kp1][2]*Aim1kp1;
                        nbim1kp1 = -1;
                    }
                    if (i==ni-1)
                    {
                        linSystemX.b_[nbijk] -= vertices_[ip1][0]*Aip1;
                        linSystemY.b_[nbijk] -= vertices_[ip1][1]*Aip1;
                        linSystemZ.b_[nbijk] -= vertices_[ip1][2]*Aip1;
                        nbip1 = -1;
                        linSystemX.b_[nbijk] -= vertices_[ip1jm1][0]*Aip1jm1;
                        linSystemY.b_[nbijk] -= vertices_[ip1jm1][1]*Aip1jm1;
                        linSystemZ.b_[nbijk] -= vertices_[ip1jm1][2]*Aip1jm1;
                        nbip1jm1 = -1;
                        linSystemX.b_[nbijk] -= vertices_[ip1jp1][0]*Aip1jp1;
                        linSystemY.b_[nbijk] -= vertices_[ip1jp1][1]*Aip1jp1;
                        linSystemZ.b_[nbijk] -= vertices_[ip1jp1][2]*Aip1jp1;
                        nbip1jp1 = -1;
                        linSystemX.b_[nbijk] -= vertices_[ip1km1][0]*Aip1km1;
                        linSystemY.b_[nbijk] -= vertices_[ip1km1][1]*Aip1km1;
                        linSystemZ.b_[nbijk] -= vertices_[ip1km1][2]*Aip1km1;
                        nbip1km1 = -1;
                        linSystemX.b_[nbijk] -= vertices_[ip1kp1][0]*Aip1kp1;
                        linSystemY.b_[nbijk] -= vertices_[ip1kp1][1]*Aip1kp1;
                        linSystemZ.b_[nbijk] -= vertices_[ip1kp1][2]*Aip1kp1;
                        nbip1kp1 = -1;
                    }
                    if (j==1)
                    {
                        linSystemX.b_[nbijk] -= vertices_[jm1][0]*Ajm1;
                        linSystemY.b_[nbijk] -= vertices_[jm1][1]*Ajm1;
                        linSystemZ.b_[nbijk] -= vertices_[jm1][2]*Ajm1;
                        nbjm1 = -1;
                        if(nbim1jm1>=0)
                        {
                            linSystemX.b_[nbijk] -= vertices_[im1jm1][0]*Aim1jm1;
                            linSystemY.b_[nbijk] -= vertices_[im1jm1][1]*Aim1jm1;
                            linSystemZ.b_[nbijk] -= vertices_[im1jm1][2]*Aim1jm1;
                            nbim1jm1 = -1;
                        }
                        if(nbip1jm1>=0)
                        {
                            linSystemX.b_[nbijk] -= vertices_[ip1jm1][0]*Aip1jm1;
                            linSystemY.b_[nbijk] -= vertices_[ip1jm1][1]*Aip1jm1;
                            linSystemZ.b_[nbijk] -= vertices_[ip1jm1][2]*Aip1jm1;
                            nbip1jm1 = -1;
                        }
                        linSystemX.b_[nbijk] -= vertices_[jm1km1][0]*Ajm1km1;
                        linSystemY.b_[nbijk] -= vertices_[jm1km1][1]*Ajm1km1;
                        linSystemZ.b_[nbijk] -= vertices_[jm1km1][2]*Ajm1km1;
                        nbjm1km1 = -1;
                        
                        linSystemX.b_[nbijk] -= vertices_[jm1kp1][0]*Ajm1kp1;
                        linSystemY.b_[nbijk] -= vertices_[jm1kp1][1]*Ajm1kp1;
                        linSystemZ.b_[nbijk] -= vertices_[jm1kp1][2]*Ajm1kp1;
                        nbjm1kp1 = -1;
                    }
                    if (j==nj-1)
                    {
                        linSystemX.b_[nbijk] -= vertices_[jp1][0]*Ajp1;
                        linSystemY.b_[nbijk] -= vertices_[jp1][1]*Ajp1;
                        linSystemZ.b_[nbijk] -= vertices_[jp1][2]*Ajp1;
                        nbjp1 = -1;
                        if(nbim1jp1>=0)
                        {
                            linSystemX.b_[nbijk] -= vertices_[im1jp1][0]*Aim1jp1;
                            linSystemY.b_[nbijk] -= vertices_[im1jp1][1]*Aim1jp1;
                            linSystemZ.b_[nbijk] -= vertices_[im1jp1][2]*Aim1jp1;
                            nbim1jp1 = -1;
                        }
                        if(nbip1jp1>=0)
                        {
                            linSystemX.b_[nbijk] -= vertices_[ip1jp1][0]*Aip1jp1;
                            linSystemY.b_[nbijk] -= vertices_[ip1jp1][1]*Aip1jp1;
                            linSystemZ.b_[nbijk] -= vertices_[ip1jp1][2]*Aip1jp1;
                            nbip1jp1 = -1;
                        }
                        linSystemX.b_[nbijk] -= vertices_[jp1km1][0]*Ajp1km1;
                        linSystemY.b_[nbijk] -= vertices_[jp1km1][1]*Ajp1km1;
                        linSystemZ.b_[nbijk] -= vertices_[jp1km1][2]*Ajp1km1;
                        nbjp1km1 = -1;
                        linSystemX.b_[nbijk] -= vertices_[jp1kp1][0]*Ajp1kp1;
                        linSystemY.b_[nbijk] -= vertices_[jp1kp1][1]*Ajp1kp1;
                        linSystemZ.b_[nbijk] -= vertices_[jp1kp1][2]*Ajp1kp1;
                        nbjp1kp1 = -1;
                    }
                    if (k==1)
                    {
                        linSystemX.b_[nbijk] -= vertices_[km1][0]*Akm1;
                        linSystemY.b_[nbijk] -= vertices_[km1][1]*Akm1;
                        linSystemZ.b_[nbijk] -= vertices_[km1][2]*Akm1;
                        nbkm1 = -1;
                        if(nbim1km1>=0)
                        {
                            linSystemX.b_[nbijk] -= vertices_[im1km1][0]*Aim1km1;
                            linSystemY.b_[nbijk] -= vertices_[im1km1][1]*Aim1km1;
                            linSystemZ.b_[nbijk] -= vertices_[im1km1][2]*Aim1km1;
                            nbim1km1 = -1;
                        }
                        if(nbip1km1>=0)
                        {
                            linSystemX.b_[nbijk] -= vertices_[ip1km1][0]*Aip1km1;
                            linSystemY.b_[nbijk] -= vertices_[ip1km1][1]*Aip1km1;
                            linSystemZ.b_[nbijk] -= vertices_[ip1km1][2]*Aip1km1;
                            nbip1km1 = -1;
                        }
                        if(nbjm1km1>=0)
                        {
                            linSystemX.b_[nbijk] -= vertices_[jm1km1][0]*Ajm1km1;
                            linSystemY.b_[nbijk] -= vertices_[jm1km1][1]*Ajm1km1;
                            linSystemZ.b_[nbijk] -= vertices_[jm1km1][2]*Ajm1km1;
                            nbjm1km1 = -1;
                        }
                        if(nbjp1km1>=0)
                        {
                            linSystemX.b_[nbijk] -= vertices_[jp1km1][0]*Ajp1km1;
                            linSystemY.b_[nbijk] -= vertices_[jp1km1][1]*Ajp1km1;
                            linSystemZ.b_[nbijk] -= vertices_[jp1km1][2]*Ajp1km1;
                            nbjp1km1 = -1;
                        }
                    }
                    if (k==nk-1)
                    {
                        linSystemX.b_[nbijk] -= vertices_[kp1][0]*Akp1;
                        linSystemY.b_[nbijk] -= vertices_[kp1][1]*Akp1;
                        linSystemZ.b_[nbijk] -= vertices_[kp1][2]*Akp1;
                        nbkp1 = -1;
                        if(nbim1kp1>=0)
                        {
                            linSystemX.b_[nbijk] -= vertices_[im1kp1][0]*Aim1kp1;
                            linSystemY.b_[nbijk] -= vertices_[im1kp1][1]*Aim1kp1;
                            linSystemZ.b_[nbijk] -= vertices_[im1kp1][2]*Aim1kp1;
                            nbim1kp1 = -1;
                        }
                        if(nbip1kp1>=0)
                        {
                            linSystemX.b_[nbijk] -= vertices_[ip1kp1][0]*Aip1kp1;
                            linSystemY.b_[nbijk] -= vertices_[ip1kp1][1]*Aip1kp1;
                            linSystemZ.b_[nbijk] -= vertices_[ip1kp1][2]*Aip1kp1;
                            nbip1kp1 = -1;
                        }
                        if(nbjm1kp1>=0)
                        {
                            linSystemX.b_[nbijk] -= vertices_[jm1kp1][0]*Ajm1kp1;
                            linSystemY.b_[nbijk] -= vertices_[jm1kp1][1]*Ajm1kp1;
                            linSystemZ.b_[nbijk] -= vertices_[jm1kp1][2]*Ajm1kp1;
                            nbjm1kp1 = -1;
                        }
                        if(nbjp1kp1>=0)
                        {
                            linSystemX.b_[nbijk] -= vertices_[jp1kp1][0]*Ajp1kp1;
                            linSystemY.b_[nbijk] -= vertices_[jp1kp1][1]*Ajp1kp1;
                            linSystemZ.b_[nbijk] -= vertices_[jp1kp1][2]*Ajp1kp1;
                            nbjp1kp1 = -1;
                        }
                    }
                    if(nbjm1km1>=0) 
                    {
                        linSystemX.A_.appendRowElement(nbijk, nbjm1km1, Ajm1km1);
                    }
                    if(nbim1km1>=0) 
                    {
                        linSystemX.A_.appendRowElement(nbijk, nbim1km1, Aim1km1);
                    }
                    if(nbkm1>=0) 
                    {
                        linSystemX.A_.appendRowElement(nbijk, nbkm1, Akm1);
                    }
                    if(nbip1km1>=0) 
                    {
                        linSystemX.A_.appendRowElement(nbijk, nbip1km1, Aip1km1);
                    }
                    if(nbjp1km1>=0) 
                    {
                        linSystemX.A_.appendRowElement(nbijk, nbjp1km1, Ajp1km1);
                    }
                    if(nbim1jm1>=0) 
                    {
                        linSystemX.A_.appendRowElement(nbijk, nbim1jm1, Aim1jm1);
                    }
                    if(nbjm1>=0) 
                    {
                        linSystemX.A_.appendRowElement(nbijk, nbjm1, Ajm1);
                    }
                    if(nbip1jm1>=0) 
                    {
                        linSystemX.A_.appendRowElement(nbijk, nbip1jm1, Aip1jm1);
                    }

                    if(nbim1>=0) 
                    {
                        linSystemX.A_.appendRowElement(nbijk, nbim1, Aim1);
                    }
                    
                    linSystemX.A_.appendRowElement(nbijk, nbijk, Aijk);
                    
                    if(nbip1>=0) 
                    {
                        linSystemX.A_.appendRowElement(nbijk, nbip1, Aip1);
                    }
                    if(nbim1jp1>=0) 
                    {
                        linSystemX.A_.appendRowElement(nbijk, nbim1jp1, Aim1jp1);
                    }
                    if(nbjp1>=0) 
                    {
                        linSystemX.A_.appendRowElement(nbijk, nbjp1, Ajp1);
                    }
                    if(nbip1jp1>=0) 
                    {
                        linSystemX.A_.appendRowElement(nbijk, nbip1jp1, Aip1jp1);
                    }
                    if(nbjm1kp1>=0) 
                    {
                        linSystemX.A_.appendRowElement(nbijk, nbjm1kp1, Ajm1kp1);
                    }
                    if(nbim1kp1>=0) 
                    {
                        linSystemX.A_.appendRowElement(nbijk, nbim1kp1, Aim1kp1);
                    }
                    if(nbkp1>=0)
                    {
                        linSystemX.A_.appendRowElement(nbijk, nbkp1, Akp1);
                    }
                    if(nbip1kp1>=0)
                    {
                        linSystemX.A_.appendRowElement(nbijk, nbip1kp1, Aip1kp1);
                    }
                    if(nbjp1kp1>=0) 
                    {
                        linSystemX.A_.appendRowElement(nbijk, nbjp1kp1, Ajp1kp1);
                    }
                }
            }
        }

        // Copy X's matrix to Y and Z. The coeffs are the same
        linSystemY.A_ = linSystemX.A_;
        linSystemZ.A_ = linSystemX.A_;
        
	    // solve the linear systems
        bool solvedX = amg.solve(linSystemX);
        if (!solvedX) Info << "Could not converge X" << endl;
        bool solvedY = amg.solve(linSystemY);
        if (!solvedY) Info << "Could not converge Y" << endl;
        bool solvedZ = amg.solve(linSystemZ);
        if (!solvedZ) Info << "Could not converge Z" << endl;
    }
    
    // Grid solved and convered. Copy back to vertices_
    for (label k = 1; k < nk; k++)
    {
        for (label j = 1; j < nj; j++)
        {
            for (label i = 1; i < ni; i++)
            {
                label ijk = vtxLabel(i, j, k);
                label nbijk = vtxLabelNoBoundary(i, j, k);
                vertices_[ijk][0] = linSystemX.x_[nbijk];
                vertices_[ijk][1] = linSystemY.x_[nbijk];
                vertices_[ijk][2] = linSystemZ.x_[nbijk];
            }
        }
    }
}


void Foam::block::createCells() const
{
    const label ni = meshDensity().x();
    const label nj = meshDensity().y();
    const label nk = meshDensity().z();

    //
    // generate cells
    //
    cells_.clear();
    cells_.setSize(nCells());

    label cellNo = 0;

    for (label k = 0; k < nk; k++)
    {
        for (label j = 0; j < nj; j++)
        {
            for (label i = 0; i < ni; i++)
            {
                cells_[cellNo].setSize(8);

                cells_[cellNo][0] =  vtxLabel(i, j, k);
                cells_[cellNo][1] =  vtxLabel(i+1, j, k);
                cells_[cellNo][2] =  vtxLabel(i+1, j+1, k);
                cells_[cellNo][3] =  vtxLabel(i, j+1, k);
                cells_[cellNo][4] =  vtxLabel(i, j, k+1);
                cells_[cellNo][5] =  vtxLabel(i+1, j, k+1);
                cells_[cellNo][6] =  vtxLabel(i+1, j+1, k+1);
                cells_[cellNo][7] =  vtxLabel(i, j+1, k+1);
                cellNo++;
            }
        }
    }
}


void Foam::block::createBoundary() const
{
    const label ni = meshDensity().x();
    const label nj = meshDensity().y();
    const label nk = meshDensity().z();

    //
    // generate boundaries on each side of the hex
    //
    boundaryPatches_.clear();
    boundaryPatches_.setSize(6);


    // x-direction

    label wallLabel = 0;
    label wallCellLabel = 0;

    // x-min
    boundaryPatches_[wallLabel].setSize(nj*nk);
    for (label k = 0; k < nk; k++)
    {
        for (label j = 0; j < nj; j++)
        {
            boundaryPatches_[wallLabel][wallCellLabel].setSize(4);

            // set the points
            boundaryPatches_[wallLabel][wallCellLabel][0] =
                vtxLabel(0, j, k);
            boundaryPatches_[wallLabel][wallCellLabel][1] =
                vtxLabel(0, j, k + 1);
            boundaryPatches_[wallLabel][wallCellLabel][2] =
                vtxLabel(0, j + 1, k + 1);
            boundaryPatches_[wallLabel][wallCellLabel][3] =
                vtxLabel(0, j + 1, k);

            // update the counter
            wallCellLabel++;
        }
    }

    // x-max
    wallLabel++;
    wallCellLabel = 0;

    boundaryPatches_[wallLabel].setSize(nj*nk);

    for (label k = 0; k < nk; k++)
    {
        for (label j = 0; j < nj; j++)
        {
            boundaryPatches_[wallLabel][wallCellLabel].setSize(4);

            // set the points
            boundaryPatches_[wallLabel][wallCellLabel][0] =
                vtxLabel(ni, j, k);
            boundaryPatches_[wallLabel][wallCellLabel][1] =
                vtxLabel(ni, j+1, k);
            boundaryPatches_[wallLabel][wallCellLabel][2] =
                vtxLabel(ni, j+1, k+1);
            boundaryPatches_[wallLabel][wallCellLabel][3] =
                vtxLabel(ni, j, k+1);

            // update the counter
            wallCellLabel++;
        }
    }

    // y-direction

    // y-min
    wallLabel++;
    wallCellLabel = 0;

    boundaryPatches_[wallLabel].setSize(ni*nk);
    for (label i = 0; i < ni; i++)
    {
        for (label k = 0; k < nk; k++)
        {
            boundaryPatches_[wallLabel][wallCellLabel].setSize(4);

            // set the points
            boundaryPatches_[wallLabel][wallCellLabel][0] =
                vtxLabel(i, 0, k);
            boundaryPatches_[wallLabel][wallCellLabel][1] =
                vtxLabel(i + 1, 0, k);
            boundaryPatches_[wallLabel][wallCellLabel][2] =
                vtxLabel(i + 1, 0, k + 1);
            boundaryPatches_[wallLabel][wallCellLabel][3] =
                vtxLabel(i, 0, k + 1);

            // update the counter
            wallCellLabel++;
        }
    }

    // y-max
    wallLabel++;
    wallCellLabel = 0;

    boundaryPatches_[wallLabel].setSize(ni*nk);

    for (label i = 0; i < ni; i++)
    {
        for (label k = 0; k < nk; k++)
        {
            boundaryPatches_[wallLabel][wallCellLabel].setSize(4);

            // set the points
            boundaryPatches_[wallLabel][wallCellLabel][0] =
                vtxLabel(i, nj, k);
            boundaryPatches_[wallLabel][wallCellLabel][1] =
                vtxLabel(i, nj, k + 1);
            boundaryPatches_[wallLabel][wallCellLabel][2] =
                vtxLabel(i + 1, nj, k + 1);
            boundaryPatches_[wallLabel][wallCellLabel][3] =
                vtxLabel(i + 1, nj, k);

            // update the counter
            wallCellLabel++;
        }
    }

    // z-direction

    // z-min
    wallLabel++;
    wallCellLabel = 0;

    boundaryPatches_[wallLabel].setSize(ni*nj);

    for (label i = 0; i < ni; i++)
    {
        for (label j = 0; j < nj; j++)
        {
            boundaryPatches_[wallLabel][wallCellLabel].setSize(4);

            // set the points
            boundaryPatches_[wallLabel][wallCellLabel][0] =
                vtxLabel(i, j, 0);
            boundaryPatches_[wallLabel][wallCellLabel][1] =
                vtxLabel(i, j + 1, 0);
            boundaryPatches_[wallLabel][wallCellLabel][2] =
                vtxLabel(i + 1, j + 1, 0);
            boundaryPatches_[wallLabel][wallCellLabel][3] =
                vtxLabel(i + 1, j, 0);

            // update the counter
            wallCellLabel++;
        }
    }

    // z-max
    wallLabel++;
    wallCellLabel = 0;

    boundaryPatches_[wallLabel].setSize(ni*nj);

    for (label i = 0; i < ni; i++)
    {
        for (label j = 0; j < nj; j++)
        {
            boundaryPatches_[wallLabel][wallCellLabel].setSize(4);

            // set the points
            boundaryPatches_[wallLabel][wallCellLabel][0] =
                vtxLabel(i, j, nk);
            boundaryPatches_[wallLabel][wallCellLabel][1] =
                vtxLabel(i + 1, j, nk);
            boundaryPatches_[wallLabel][wallCellLabel][2] =
                vtxLabel(i + 1, j + 1, nk);
            boundaryPatches_[wallLabel][wallCellLabel][3] =
                vtxLabel(i, j + 1, nk);

            // update the counter
            wallCellLabel++;
        }
    }
}


void Foam::block::clearGeom()
{
    vertices_.clear();
    cells_.clear();
    boundaryPatches_.clear();
}


// ************************************************************************* //
