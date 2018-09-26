/** A template for using Peer (Group) Execution.
@author Christopher Ferrall
**/
#include "PeerCommunicationTestA.h"

MyPeer::MyPeer() {
    Peer();
    }

MyPeer::Run() {
  Volume = LOUD;
  decl me = new MyPeer();
  me.Buffer = matrix(ID+1);
  me->Sum(1);
  if (IamClient) println("Sum of node ranks + 1 = ",me.Buffer[0],".  Answer should be ",Nodes*(Nodes+1)/2);
  delete me;
  }
