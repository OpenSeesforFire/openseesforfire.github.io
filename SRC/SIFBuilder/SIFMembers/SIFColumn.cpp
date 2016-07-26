#include <SIFColumn.h>

SIFColumn::SIFColumn(int tag, int typeTag, int jt1, int jt2,const ID& theMemInfo, double gamma):SIFMember(tag, 3 ,typeTag,theMemInfo)
{
  //typeTag indicates composite beam or single beam

  //it can be extended to concrete filled steel column and other column types
  
  ConnectedJoints = new ID(2);
  ConnectedJoints->resize(2);
  (*ConnectedJoints)(0) = jt1;
  (*ConnectedJoints)(1) = jt2;
}




SIFColumn::~SIFColumn()
{
	//
}



