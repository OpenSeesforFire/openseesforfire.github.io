#include <SIFCompartment3D.h>

SIFCompartment3D::SIFCompartment3D(int tag, double cxx1, double cxx2, double cxx3, double cxx4, double rbxx1, double rbxx2, double rbxx3, double rbxx4, double rsx1, int wll, int wlr, double wlf, double wlb):TaggedObject(tag),Mem(0)
{
		Mem = new Vector(9);
		(*Mem)(0) = cxx1;
		(*Mem)(1) = cxx2;
		(*Mem)(2) = cxx3;
		(*Mem)(3) = cxx4;
		(*Mem)(4) = rbxx1;
		(*Mem)(5) = rbxx2;
		(*Mem)(6) = rbxx3;
		(*Mem)(7) = rbxx4;
		(*Mem)(8) = rsx1;
		(*Mem)(9) = wll;
		(*Mem)(10) = wlr;
		(*Mem)(11) = wlf;
		(*Mem)(12) = wlb;
		
		
}

SIFCompartment3D::~SIFCompartment3D()
{
	//
}

const ID&
SIFCompartment3D::getConnectedXBeams(void){

	ID connectedMembers =0;
	return connectedMembers;
}


const ID&
SIFCompartment3D::getConnectedYBeams(void){

	ID connectedMembers =0;
	return connectedMembers;
}

const ID&
SIFCompartment3D::getConnectedColumns(void){

	ID connectedMembers =0;
	return connectedMembers;
}

const ID&
SIFCompartment3D::getConnectedSlabs(void){

	ID connectedSlabs =0;
	return connectedSlabs;
}