#include "DrawTree.oxh"

 /** This function is used to create a game tree for a 2 x 2 matrix in a sequential game.
@param p1  n x m matrix of payoffs for player 1 (player strategies are rows)
@param p22  n x m matrix of payoffs for player 2 (player strategies are columns)
@param fn string, file name to save tree to (default=<q>Tree</q>).

The game tree will save in the directory where the package has been saved to

**/
DrawTree(p1, p2, fn){
    SetDrawWindow("Sequential Game");
    p2 = p2';
    // Create the plot area
    DrawAxis(0,1,0,0,30,0,0,0,0);
    DrawAxis(0,0,0,0,20,0,0,0,0);
    //Draw all the lines
    DrawLine(0, 15,20,10,15,1); //Top left
    DrawLine(0, 15,20,20,15,1);	//Top Right

    DrawLine(0, 10,14,7,9,1);  //AA
    DrawLine(0, 10,14,13,9,1); //AB

    DrawLine(0, 20,14,17,9,1); //BA
    DrawLine(0, 20,14,23,9,1); //BB

    //Place text
    DrawText(0, "1", 14.81,20.25);	  //Top
    DrawText(0, "2", 9.81,14.25);	  //Left
    DrawText(0, "2", 19.81,14.25);	  //Right
    DrawText(0, "A", 11.5,17.75); 	  //Top A
    DrawText(0, "B", 18,17.75);	      //Top B
    DrawText(0, "D", 7.75,11.75);	  //AA
    DrawText(0, "E", 11.8,11.65);	  //AB
    DrawText(0, "D", 17.75,11.75);	  //BA
    DrawText(0, "E", 21.8,11.65);	  //BB

    //Place Circles
    DrawSymbol(0, 14.5, 20, 15.5, 21, PL_CIRCLE, 1); //Top
    DrawSymbol(0, 9.5, 14, 10.5, 15, PL_CIRCLE, 1);	 //Left
    DrawSymbol(0, 19.5, 14, 20.5, 15, PL_CIRCLE, 1); //Right

    //Place Payoffs
    //decl p1 = <1, 2; 3, 4>, p2 = <5,6;7,8>;

    DrawText(0, sprint("( ",p1[0][0]," , ",p2[0][0]," )"),6,8.25);  //AA
    DrawText(0, sprint("( ",p1[0][1]," , ",p2[0][1]," )"),12,8.25); //AB
    DrawText(0, sprint("( ",p1[1][0]," , ",p2[1][0]," )"),16,8.25); //BA
    DrawText(0, sprint("( ",p1[1][1]," , ",p2[1][1]," )"),22,8.25); //BB

    SaveDrawWindow(fn+".pdf");
}
