        newPs = (*rmesh).points();
//- set-up interpolator
	volPointInterpolation vpi(*rmesh);

//- perform interpolation
	pointScalarField deltaIcep = vpi.interpolate(filmSWIM->deltaIce()-olddeltaIce);
//	Info<<"new value deltaIce:"<<endl;
//	Info<<filmSWIM->deltaIce()<<endl;
//	Info<<"old value deltaIce:"<<endl;
//	Info<<deltaIcep<<endl;
	const labelList& intfilmPoints=(*rmesh).boundaryMesh()[intfilmPatchID].meshPoints();
	const labelListList& intfilmpointFaces=(*rmesh).boundaryMesh()[intfilmPatchID].pointFaces();
	const vectorField& intfilmnf=(*rmesh).boundary()[intfilmPatchID].Sf();

	forAll(intfilmPoints, index) 
	{
	    label i=intfilmPoints[index];
	    if (deltaIcep[i]>0.0) {
//	    Info<<"index="<<index<<" i="<<i<<" deltaIce="<<deltaIcep[i]<<endl;
		vector d(0.0, 0.0, 0.0);
		forAll(intfilmpointFaces[index], faceI) {
//		  Info<<"faceI="<<faceI<<" nf="<<intfilmnf[intfilmpointFaces[index][faceI]]<<" "<<intfilmpointFaces[index][faceI]<<" <="<<intfilmnf.size()<<endl;
		  d+=intfilmnf[intfilmpointFaces[index][faceI]]/mag(intfilmnf[intfilmpointFaces[index][faceI]]);
		}
//		Info<<"d="<<d<<endl;
//		Info<<"Point="<<newPs[i]<<endl;
//		newPs[i]+=deltaIcep[i]*d/mag(d);
		newPs[i].x()+=deltaIcep[i]*d.x()/mag(d);
		newPs[i].z()+=deltaIcep[i]*d.z()/mag(d);
    	    }
	}
	
	const labelList& inttopfilmPoints=(*rmesh).boundaryMesh()[inttopfilmPatchID].meshPoints();
	const labelListList& inttopfilmpointFaces=(*rmesh).boundaryMesh()[inttopfilmPatchID].pointFaces();
	const vectorField& inttopfilmnf=(*rmesh).boundary()[inttopfilmPatchID].Sf();
	forAll(inttopfilmPoints, index) 
	{
	    label i=inttopfilmPoints[index];
	    if (deltaIcep[i]>0.0) {
		vector d(0.0, 0.0, 0.0);
		forAll(inttopfilmpointFaces[index], faceI) {
		  d+=inttopfilmnf[inttopfilmpointFaces[index][faceI]]/mag(inttopfilmnf[inttopfilmpointFaces[index][faceI]]);
		}
//		newPs[i]-=deltaIcep[i]*d/mag(d);
		newPs[i].x()-=deltaIcep[i]*d.x()/mag(d);
		newPs[i].z()-=deltaIcep[i]*d.z()/mag(d);
    	    }
	}

	(*rmesh).movePoints(newPs);

//- set-up interpolator
	primitivePatchInterpolation priInterpolator(rmesh->boundaryMesh()[intfilmPatchID]);
	scalarField dIcePatch = filmSWIM->deltaIce().boundaryField()[intfilmPatchID]-olddeltaIce.boundaryField()[intfilmPatchID];
	scalarField dIcePoint = priInterpolator.faceToPointInterpolate(dIcePatch);

	patchToPatchInterpolation patchInterpolator 
	(
	    rmesh->boundaryMesh()[intfilmPatchID], 
	    mesh.boundaryMesh()[filmPatchID],
	    intersection::FULL_RAY,
	    intersection::CONTACT_SPHERE
	);

//- perform interpolation
	scalarField iceValues = patchInterpolator.pointInterpolate(dIcePoint);

	mesh.moving(true);

	const labelList& filmPoints=mesh.boundaryMesh()[filmPatchID].meshPoints();
	const labelListList& filmpointFaces=mesh.boundaryMesh()[filmPatchID].pointFaces();
	const vectorField& filmnf=mesh.boundary()[filmPatchID].Sf();

	forAll(filmPoints, index) 
	{
//	    label i=filmPoints[index];
	    if (iceValues[index]>0.0) {
		vector d(0.0, 0.0, 0.0);
		forAll(filmpointFaces[index], faceI) {

		  d+=filmnf[filmpointFaces[index][faceI]]/mag(filmnf[filmpointFaces[index][faceI]]);
		}
//		pDisp[index]+=iceValues[index]*d/mag(d);
		pDisp[index].y()=0.0; 
		pDisp[index].x()-=iceValues[index]*d.x()/mag(d);
		pDisp[index].z()-=iceValues[index]*d.z()/mag(d);
		
    	    }
/*	    scalar oldPx=mesh.points()[i].x();
	    scalar oldPz=mesh.points()[i].z();
	    scalar r2=Foam::sqrt(oldPx*oldPx+oldPz*oldPz);
	    pDisp[index].y()=0.0; 
	    pDisp[index].x()+=oldPx*((iceValues[index]+r)/r2-1.0);
	    pDisp[index].z()+=oldPz*((iceValues[index]+r)/r2-1.0);
*/
	}

    	// Do any mesh changes
        parcels.storeGlobalPositions();
    	mesh.update();
    	parcels.updateMesh();
        mesh.moving(false);
