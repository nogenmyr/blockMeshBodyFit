/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | Unsupported Contributions for OpenFOAM
 \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
   \\/     M anipulation  |
-------------------------------------------------------------------------------
2013-06-01 Karl-Johan Nogenmyr. 

Constructur got longer to initialize new data members

Added implementation of:
        void calculatePointNormals();
        void projectBlockPoints();
        void createProjectedEdges();


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

#include "blockMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::blockMesh::blockMesh::verboseOutput(false);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockMesh::blockMesh(const IOdictionary& dict, const word& regionName)
:
    blockPointField_(dict.lookup("vertices")),
    scaleFactor_(1.0),
    topologyPtr_(createTopology(dict, regionName)),
    snapSurfaces_(0),
    fixedPoint_(blockPointField_.size(), true),
    projectionList_(topologyPtr_->nFaces(), -1),
    faceNormal_(topologyPtr_->nFaces(), Foam::vector(0,0,0)),
    blockPointProjectionList_(topologyPtr_->nPoints(), -1),
    pointNormals_(topologyPtr_->nPoints(), Foam::vector(0,0,0)),
    searchLength_(-1)
{
    
    if (dict.found("snapFaces"))
    {
        Info<< nl << "Reading snapFaces section" << endl;
        searchLength_ = readScalar(dict.lookup("searchLength"));
        readSnapFaces(dict);
        calculatePointNormals();
        forAll(edges_, edgeI) //pre defined edges should not move...
        {
            curvedEdge& ed = edges_[edgeI];
            fixedPoint_[ed.start()] = true;
            fixedPoint_[ed.end()] = true;
        }
        projectBlockPoints();
        pointNormals_ *= 0;
        calculatePointNormals(); //Call again after projectBlockPoints()
        createProjectedEdges();
    }
    
    blockList& blocks = *this;
    
    // after moving block's corners, we have to update edges also
    forAll(blocks, blockI) 
    {
        blocks[blockI].makeBlockEdges();
    }
    
    calcMergeInfo();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blockMesh::~blockMesh()
{
    delete topologyPtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockMesh::calculatePointNormals()
{
    vectorField centreNormals(topologyPtr_->nFaces(), Foam::vector(0,0,0));

    forAll(topologyPtr_->points(), pointI)
    {
        forAll(projectionList_, faceI)
        {
            if(projectionList_[faceI] >= 0)
            {
                forAll(topologyPtr_->faces()[faceI], fpI)
                {
                    if(topologyPtr_->faces()[faceI][fpI] == pointI)
                    {
                        scalar dist = mag
                        (
                            blockPointField_[pointI]
                          - topologyPtr_->faces()[faceI].centre(blockPointField_)
                        );
                         
                        Foam::vector faceNormal = 
                         (
                            topologyPtr_->faces()[faceI].normal(blockPointField_)
                        );

                        if((faceNormal & faceNormal_[faceI]) < 0.) faceNormal *= -1;
                        pointNormals_[pointI] += faceNormal/(dist*mag(faceNormal));

                        blockPointProjectionList_[pointI] = projectionList_[faceI];
                    }
                }
            }
        }
                    
        if(blockPointProjectionList_[pointI]>=0)
        {
            pointNormals_[pointI] /= mag(pointNormals_[pointI]);
            fixedPoint_[pointI] = false;
        }
    }
}

void Foam::blockMesh::projectBlockPoints()
{
    pointField newPoints = blockPointField_;

    forAll(topologyPtr_->points(), pointI)
    {
        label surfI = blockPointProjectionList_[pointI];
        if(surfI >= 0 && !fixedPoint_[pointI])
        {
            point& p = newPoints[pointI];
            pointField start(2, p);
            pointField end(2);
            end[0] = p+searchLength_*pointNormals_[pointI];
            end[1] = p-searchLength_*pointNormals_[pointI];
            List<pointIndexHit> hits(2);
            snapSurfaces_[surfI].findLine(start, end, hits);
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
    }
    blockPointField_ = newPoints;
    topologyPtr_->movePoints(newPoints);
}

void Foam::blockMesh::createProjectedEdges()
{
    label nEdges = edges_.size();
    forAll(topologyPtr_->edges(), edgeI)
    {
        const edge& ed = topologyPtr_->edges()[edgeI];

        label nFacesWithEdge = 0;
        forAll(projectionList_, faceI)
        {
            if(projectionList_[faceI] >= 0)
            { 
                bool point1Found = false;
                bool point2Found = false;
                forAll(topologyPtr_->faces()[faceI], fpI)
                {
                    const label& facep = topologyPtr_->faces()[faceI][fpI];
                    if (facep == ed.start()) point1Found = true;
                    if (facep == ed.end()) point2Found = true;
                }
                if (point1Found && point2Found) nFacesWithEdge++;
            }
        }
        bool isPreDefinedEdge = false;
        forAll(edges_, curvedEdgeI)
        {
            curvedEdge& cEd = edges_[curvedEdgeI];
            if 
            (
                (cEd.start() == ed.start() || cEd.start() == ed.end())
             && (cEd.end() == ed.start() || cEd.end() == ed.end())
            )  isPreDefinedEdge = true;
        }
        if(nFacesWithEdge && !isPreDefinedEdge)
        {
            vector& normalStart(pointNormals_[ed.start()]);
            vector& normalEnd(pointNormals_[ed.end()]);
            edges_.setSize(nEdges + 1);

            edges_.set
            (
                nEdges,
                new projectionEdge
                (
                    blockPointField_, 
                    ed.start(), 
                    ed.end(), 
                    normalStart, 
                    normalEnd, 
                    snapSurfaces_,
                    searchLength_,
                    blockPointProjectionList_[ed.start()]                    
                )
            );
            nEdges++;
        }
    }
}

void Foam::blockMesh::verbose(const bool on)
{
    verboseOutput = on;
}


const Foam::pointField& Foam::blockMesh::blockPointField() const
{
    return blockPointField_;
}


const Foam::polyMesh& Foam::blockMesh::topology() const
{
    if (!topologyPtr_)
    {
        FatalErrorIn("blockMesh::topology() const")
            << "topologyPtr_ not allocated"
            << exit(FatalError);
    }

    return *topologyPtr_;
}


Foam::PtrList<Foam::dictionary> Foam::blockMesh::patchDicts() const
{
    const polyPatchList& patchTopologies = topology().boundaryMesh();

    PtrList<dictionary> patchDicts(patchTopologies.size());

    forAll(patchTopologies, patchI)
    {
        OStringStream os;
        patchTopologies[patchI].write(os);
        IStringStream is(os.str());
        patchDicts.set(patchI, new dictionary(is));
    }
    return patchDicts;
}


Foam::scalar Foam::blockMesh::scaleFactor() const
{
    return scaleFactor_;
}


const Foam::pointField& Foam::blockMesh::points() const
{
    if (points_.empty())
    {
        createPoints();
    }

    return points_;
}


const Foam::cellShapeList& Foam::blockMesh::cells() const
{
    if (cells_.empty())
    {
        createCells();
    }

    return cells_;
}


const Foam::faceListList& Foam::blockMesh::patches() const
{
    if (patches_.empty())
    {
        createPatches();
    }

    return patches_;
}


Foam::wordList Foam::blockMesh::patchNames() const
{
    return topology().boundaryMesh().names();
}


//Foam::wordList Foam::blockMesh::patchTypes() const
//{
//    return topology().boundaryMesh().types();
//}
//
//
//Foam::wordList Foam::blockMesh::patchPhysicalTypes() const
//{
//    return topology().boundaryMesh().physicalTypes();
//}


Foam::label Foam::blockMesh::numZonedBlocks() const
{
    label num = 0;

    forAll(*this, blockI)
    {
        if (operator[](blockI).zoneName().size())
        {
            num++;
        }
    }

    return num;
}


void Foam::blockMesh::writeTopology(Ostream& os) const
{
    const pointField& pts = topology().points();

    forAll(pts, pI)
    {
        const point& pt = pts[pI];

        os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
    }

    const edgeList& edges = topology().edges();

    forAll(edges, eI)
    {
        const edge& e = edges[eI];

        os << "l " << e.start() + 1 << ' ' << e.end() + 1 << endl;
    }
}

// ************************************************************************* //
