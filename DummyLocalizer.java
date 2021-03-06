package model;

import control.EstimatorInterface;
import java.util.ArrayList;

public class DummyLocalizer implements EstimatorInterface {
		
	private int rows, cols, head;
	private double [][] T;
	private ArrayList<double[]> O_diags;
    private double[] Onothing;
    private int trueRow, trueCol, trueHead;
    private int percRow, percCol;
    private double[] F;
    private double avgDist = 0;
    private int its = 0;
    private double totDist = 0;

	public DummyLocalizer( int rows, int cols, int head) {
        int n = rows*cols*head;
        this.rows = rows;
        this.cols = cols;
        this.head = head;
        double x = Math.random();
        x = x*rows;
        double y = Math.random();
        y = y*cols;
        double h = Math.random();
        h = h*head;
        trueRow = (int)x;
        trueCol = (int)y;
        trueHead = (int)h;
        T = new double[n][n];
        setCorners();
        setEdges();
        setInner();
        setONothing();
        O_diags = new ArrayList<>();
        setDiags();
        percCol=trueRow;
        percRow=trueCol;

        F = O_diags.get(trueRow*rows+trueCol);

	}
    private void setDiags(){
	    for(int i = 0;i<rows*cols;i++){
	        double[] O = new double[rows*cols*head];
	        for(int j = 0;j<rows*cols;j++){
	            int ir = i/rows;
	            int ic = i%cols;
                int jr = j/rows;
                int jc = j%cols;
                double dist = Math.sqrt(Math.abs(Math.pow(ir-jr,2)+Math.pow(ic-jc,2)));
	            if(dist<0.8){
                    O[j*4] = 0.1;
                    O[j*4+1] = 0.1;
                    O[j*4+2] = 0.1;
                    O[j*4+3] = 0.1;
                }else if(dist >0.8 && dist <1.5){
                    O[j*4] = 0.05;
                    O[j*4+1] = 0.05;
                    O[j*4+2] = 0.05;
                    O[j*4+3] = 0.05;

                }else if(dist> 1.5 && dist < 2.83){
                    O[j*4] = 0.025;
                    O[j*4+1] = 0.025;
                    O[j*4+2] = 0.025;
                    O[j*4+3] = 0.025;
                } else{
                    O[j*4] = 0;
                    O[j*4+1] = 0;
                    O[j*4+2] = 0;
                    O[j*4+3] = 0;
                }
            }
            O_diags.add(O);
        }

    }
	private void setCorners(){
        //S[0][0] ////Going east, and south
        T[0][head+1] = 0.5; //00N
        T[0][2+rows*head] = 0.5;
        T[1][head+1] = 0.7; //OOE
        T[1][2+rows*head] = 0.3;
        T[2][head+1] = 0.3; //00S
        T[2][2+rows*head] = 0.7;
        T[3][head+1] = 0.5; //00W
        T[3][2+rows*head] = 0.5;
        //S[0][rows-1] //// going west, and south
        T[head*(cols-1)][head*(cols-1)-head+3] = 0.5; //0CN
        T[head*(cols-1)][head*(cols-1)+cols*head+2] = 0.5;
        T[head*(cols-1)+1][head*(cols-1)-head+3] = 0.5; //OCE
        T[head*(cols-1)+1][head*(cols-1)+cols*head+2] = 0.5;
        T[head*(cols-1)+2][head*(cols-1)-head+3] = 0.3; //0CS
        T[head*(cols-1)+2][head*(cols-1)+cols*head+2] = 0.7;
        T[head*(cols-1)+3][head*(cols-1)-head+3] = 0.7; //0CW
        T[head*(cols-1)+3][head*(cols-1)+cols*head+2] = 0.3;

        //S[rows-1][0] ///going north, and east
        T[(rows-1)*rows*head][(rows-2)*rows*head] = 0.7; //0CN
        T[(rows-1)*rows*head][(rows-1)*rows*head+head+1] = 0.3;
        T[(rows-1)*rows*head+1][(rows-2)*rows*head] = 0.3; //OCE
        T[(rows-1)*rows*head+1][(rows-1)*rows*head+head+1] = 0.7;
        T[(rows-1)*rows*head+2][(rows-2)*rows*head] = 0.5; //0CS
        T[(rows-1)*rows*head+2][(rows-1)*rows*head+head+1] = 0.5;
        T[(rows-1)*rows*head+3][(rows-2)*rows*head] = 0.5; //0CW
        T[(rows-1)*rows*head+3][(rows-1)*rows*head+head+1] = 0.5;

        //S[rows-1][0] //going north, and west

        T[cols*head*(rows-1)+head*(cols-1)][cols*head*(rows-1)+head*(cols-1)-cols*head] = 0.7; //0CN
        T[cols*head*(rows-1)+head*(cols-1)][cols*head*(rows-1)+head*(cols-1)-head+3] = 0.3;
        T[cols*head*(rows-1)+head*(cols-1)+1][cols*head*(rows-1)+head*(cols-1)-cols*head] = 0.5; //OCE
        T[cols*head*(rows-1)+head*(cols-1)+1][cols*head*(rows-1)+head*(cols-1)-head+3] = 0.5;
        T[cols*head*(rows-1)+head*(cols-1)+2][cols*head*(rows-1)+head*(cols-1)-cols*head] = 0.5; //0CS
        T[cols*head*(rows-1)+head*(cols-1)+2][cols*head*(rows-1)+head*(cols-1)-head+3] = 0.5;
        T[cols*head*(rows-1)+head*(cols-1)+3][cols*head*(rows-1)+head*(cols-1)-cols*head] = 0.3; //0CW
        T[cols*head*(rows-1)+head*(cols-1)+3][cols*head*(rows-1)+head*(cols-1)-head+3] = 0.7;
    }
    private void setEdges(){
	    for(int i = 1;i<cols-1;i++){ //north edge
	        for(int k = 0;k<head;k++) {
                if ((i * head + k) % 4 == 0) { //N
                    T[i * head + k][(i-1) * head + 3] = 0.333;
                    T[i * head + k][(i+1) * head + 1] = 0.333;
                    T[i * head + k][(i) * head + 2 + rows*head] = 0.333;
                } else if ((i * head + k) % 4 == 1) { //e
                    T[i * head + k][(i-1) * head + 3] = 0.15;
                    T[i * head + k][(i+1) * head + 1] = 0.7;
                    T[i * head + k][(i) * head + 2 + rows*head] = 0.15;
                }else if ((i * head + k) % 4 == 2) { //s
                    T[i * head + k][(i-1) * head + 3] = 0.15;
                    T[i * head + k][(i+1) * head + 1] = 0.15;
                    T[i * head + k][(i) * head + 2 + rows*head] = 0.7;
                }else if ((i * head + k) % 4 == 3) { //w
                    T[i * head + k][(i-1) * head + 3] = 0.7;
                    T[i * head + k][(i+1) * head + 1] = 0.15;
                    T[i * head + k][(i) * head + 2 + rows*head] = 0.15;
                }
            }
        }
        for(int i = 1;i<cols-1;i++){ //east edge
            for(int k = 0;k<head;k++) {
                if ((i * rows*head + (cols-1)*head+ k) % 4 == 0) { //N
                    T[i * rows*head + (cols-1)*head+k][i * rows*head + (cols-2)*head + 3] = 0.15;//w
                    T[i * rows*head + (cols-1)*head+k][(i) * rows*head + (cols-1)*head - rows*head] = 0.7;//n
                    T[i * rows*head + (cols-1)*head+k][(i) * rows*head + (cols-1)*head + 2 + rows*head] = 0.15;//s
                } else if ((i * rows*head + (cols-1)*head+k) % 4 == 1) { //e
                    T[i * rows*head + (cols-1)*head+k][i * rows*head + (cols-2)*head + 3] = 0.333;
                    T[i * rows*head + (cols-1)*head+k][(i) * rows*head + (cols-1)*head - rows*head] = 0.333;
                    T[i * rows*head + (cols-1)*head+k][(i) * rows*head + (cols-1)*head + 2 + rows*head] = 0.333;
                }else if ((i * rows*head + (cols-1)*head+k) % 4 == 2) { //s
                    T[i * rows*head + (cols-1)*head+k][i * rows*head + (cols-2)*head + 3] = 0.15;
                    T[i * rows*head + (cols-1)*head+k][(i) * rows*head + (cols-1)*head - rows*head] = 0.15;
                    T[i * rows*head + (cols-1)*head+k][(i) * rows*head + (cols-1)*head + 2 + rows*head] = 0.7;
                }else if ((i * rows*head + (cols-1)*head+k) % 4 == 3) { //w
                    T[i * rows*head + (cols-1)*head+k][i * rows*head + (cols-2)*head + 3] = 0.7;
                    T[i * rows*head + (cols-1)*head+k][(i) * rows*head + (cols-1)*head - rows*head] = 0.15;
                    T[i * rows*head + (cols-1)*head+k][(i) * rows*head + (cols-1)*head + 2 + rows*head] = 0.15;
                }
            }
        }
        for(int i = 1;i<cols-1;i++){ //south edge
            for(int k = 0;k<head;k++) {
                if ((rows*head*(rows-1) + i * head + k) % 4 == 0) { //N
                    T[rows*head*(rows-1) + i * head + k][rows*head*(rows-1) + (i-1) * head + 3] = 0.15;
                    T[rows*head*(rows-1) + i * head + k][rows*head*(rows-1) + (i+1) * head + 1] = 0.15;
                    T[rows*head*(rows-1) + i * head + k][rows*head*(rows-1) + (i) * head - rows*head] = 0.7;
                } else if ((rows*head*(rows-1) + i * head + k) % 4 == 1) { //e
                    T[rows*head*(rows-1) + i * head + k][rows*head*(rows-1) + (i-1) * head + 3] = 0.15;
                    T[rows*head*(rows-1) + i * head + k][rows*head*(rows-1) + (i+1) * head + 1] = 0.7;
                    T[rows*head*(rows-1) + i * head + k][rows*head*(rows-1) + (i) * head - rows*head] = 0.15;
                }else if ((rows*head*(rows-1) + i * head + k) % 4 == 2) { //s
                    T[rows*head*(rows-1) + i * head + k][rows*head*(rows-1) + (i-1) * head + 3] = 0.333;
                    T[rows*head*(rows-1) + i * head + k][rows*head*(rows-1) + (i+1) * head + 1] = 0.333;
                    T[rows*head*(rows-1) + i * head + k][rows*head*(rows-1) + (i) * head - rows*head] = 0.333;
                }else if ((rows*head*(rows-1) + i * head + k) % 4 == 3) { //w
                    T[rows*head*(rows-1) + i * head + k][rows*head*(rows-1) + (i-1) * head + 3] = 0.7;
                    T[rows*head*(rows-1) + i * head + k][rows*head*(rows-1) + (i+1) * head + 1] = 0.15;
                    T[rows*head*(rows-1) + i * head + k][rows*head*(rows-1) + (i) * head - rows*head] = 0.15;
                }
            }
        }
        for(int i = 1;i<cols-1;i++){ //west edge
            for(int k = 0;k<head;k++) {
                if ((i * rows*head + k) % 4 == 0) { //N
                    T[i * rows*head+k][(i) * rows*head + cols + 1] = 0.15;//e
                    T[i * rows*head+k][(i) * rows*head - rows*head] = 0.7;//n
                    T[i * rows*head+k][(i) * rows*head + 2 + rows*head] = 0.15;//s
                } else if ((i * rows*head+k) % 4 == 1) { //e
                    T[i * rows*head+k][(i) * rows*head + cols + 1] = 0.7;
                    T[i * rows*head+k][(i) * rows*head - rows*head] = 0.15;
                    T[i * rows*head+k][(i) * rows*head  + 2 + rows*head] = 0.15;
                }else if ((i * rows*head+k) % 4 == 2) { //s
                    T[i * rows*head+k][(i) * rows*head +cols + 1] = 0.15;
                    T[i * rows*head+k][(i) * rows*head- rows*head] = 0.15;
                    T[i * rows*head+k][(i) * rows*head+ 2 + rows*head] = 0.7;
                }else if ((i * rows*head+k) % 4 == 3) { //w
                    T[i * rows*head+k][(i) * rows*head +cols+ 1] = 0.333;
                    T[i * rows*head+k][(i) * rows*head- rows*head] = 0.333;
                    T[i * rows*head+k][(i) * rows*head + 2 + rows*head] = 0.333;
                }
            }
        }

    }
    private void setInner(){
        for (int i = 1;i<rows-1;i++){
            for (int j = 1;j<cols-1;j++){
                for(int k = 0;k<head;k++){
                    if((i*rows*head+j*head+k)%4==0) { //N
                        T[i * rows * head + j * head+k][(i - 1) * rows * head + j * head] = 0.7;
                        T[i * rows * head + j * head+k][i * rows * head + (j - 1) * head+3] = 0.1;
                        T[i * rows * head + j * head+k][(i + 1) * rows * head + j * head+2] = 0.1;
                        T[i * rows * head + j * head+k][i * rows * head + (j + 1) * head+1] = 0.1;
                    } else if((i*rows*head+j*head+k)%4==1) { //E
                        T[i * rows * head + j * head+k][(i - 1) * rows * head + j * head] = 0.1;
                        T[i * rows * head + j * head+k][i * rows * head + (j - 1) * head+3] = 0.1;
                        T[i * rows * head + j * head+k][(i + 1) * rows * head + j * head+2] = 0.1;
                        T[i * rows * head + j * head+k][i * rows * head + (j + 1) * head+1] = 0.7;
                    }else if((i*rows*head+j*head+k)%4==2) { //S
                        T[i * rows * head + j * head+k][(i - 1) * rows * head + j * head] = 0.1;
                        T[i * rows * head + j * head+k][i * rows * head + (j - 1) * head+3] = 0.1;
                        T[i * rows * head + j * head+k][(i + 1) * rows * head + j * head+2] = 0.7;
                        T[i * rows * head + j * head+k][i * rows * head + (j + 1) * head+1] = 0.1;
                    }else if((i*rows*head+j*head+k)%4==3) { //W
                        T[i * rows * head + j * head + k][(i - 1) * rows * head + j * head] = 0.1;
                        T[i * rows * head + j * head + k][i * rows * head + (j - 1) * head + 3] = 0.7;
                        T[i * rows * head + j * head + k][(i + 1) * rows * head + j * head + 2] = 0.1;
                        T[i * rows * head + j * head + k][i * rows * head + (j + 1) * head + 1] = 0.1;
                    }
                }
            }
        }
    }
	public int getNumRows() {
		return rows;
	}
	public int getNumCols() {
		return cols;
	}
	public int getNumHead() {
		return head;
	}
	public double getTProb( int x, int y, int h, int nX, int nY, int nH) {
        return T[cols*head*x+head*y+h][cols*head*nX + head*nY + nH];
	}
	public double getOrXY( int rX, int rY, int x, int y, int h) {

	    if(rX == -1){
	        return Onothing[head*rows*x+y*head+h];
        }
		return O_diags.get(rows*rX + rY)[head*rows*x+y*head+h];
	}
	public int[] getCurrentTrueState() {
		int[] ret = new int[3];
		ret[0] = trueRow;
		ret[1] = trueCol;
		ret[2] = trueHead;
		return ret;

	}
	public int[] getCurrentReading() {
		int[] ret = new int[2];
		ret[0] = percRow;
		ret[1] = percCol;
		return ret;
	}
	public double getCurrentProb( int x, int y) {
		double sum = F[rows*x*head+head*y] + F[rows*x*head+head*y+1] +  F[rows*x*head+head*y+2] + F[rows*x*head+head*y+3];
		return sum;////O_diags.get(rows*x+y)[head*rows*trueRow+head*trueCol+trueHead];
	}
	public void update() {
        // gå ett steg enl T
        takeStep();
        // uppdatera perceived location efter fördelningen
        updateReading();
        // räkna ut nya F
        double[] O_diag;
        if(percRow == -1){
            O_diag = Onothing;
        } else{
            O_diag = O_diags.get(percRow*rows+percCol);
        }

        double[][] O_mat = new double[rows*cols*head][rows*cols*head];
        for(int i = 0;i<rows*cols*head;i++){
            for(int j = 0;j<rows*cols*head;j++){
                if(i==j){
                    O_mat[i][j] = O_diag[i];
                }
            }
        }
        F = forwardAlgorithm(F, O_mat);
        
        
        //räkna ut manhattan distance och printa iteration och snitt
        its++;
        double[] fProbs = new double[rows*cols];
        for(int i = 0;i<rows*cols;i++) {
        	fProbs[i] = getCurrentProb(i/rows, i%cols);
        }
        int maxIndex = 0;
        double max = -1;
        for(int i = 0;i< fProbs.length;i++) {
        	max = Math.max(max, fProbs[i]);
        	if(max == fProbs[i]) {
        		maxIndex = i;
        	}
        }
        int maxRow = maxIndex/rows;
        int maxCol = maxIndex%cols;
        int manHDist = Math.abs(maxRow-trueRow) + Math.abs(maxCol - trueCol);
        totDist += manHDist;
        avgDist = totDist/its;
        System.out.print("iteration: " + its);
        System.out.println(", average distance: "+avgDist);
        
	}
	private void updateReading(){
		int pos = trueRow*rows+trueCol;
        percRow = -1;
        percCol = -1;
        double[] d = new double[O_diags.get(pos).length/4+1];
        double[] od = O_diags.get(pos);
        for(int i = 0;i<d.length-1;i++) {
        	d[i] = od[4*i]*4;
        }
        d[d.length-1] = Onothing[pos];
        for(double nd :d){
            //System.out.print(nd + "  ");
        }
        int ans;
        double a = Math.random();
        double acc = 0;
        int i = 0;
        while(true){
            if(a<d[i]+acc){
                ans = i;
                break;
            }else{
                acc+=d[i];
                i++;
            }
        }
        percRow = ans/rows;
        percCol = ans%cols;
    }
    private double[] forwardAlgorithm(double[] F, double[][] O) {
        double [][] tempF = new double[1][F.length];
        tempF[0] = F;
        double[][] FTranspose = transposeMatrix(tempF);
        double[][] TTranspose = transposeMatrix(T);

        double[][] FMatrix = transposeMatrix(mulMatrices(mulMatrices(O,TTranspose), FTranspose));
        double sum = 0;
        for(int i = 0;i< FMatrix[0].length;i++){
            sum +=FMatrix[0][i];
        }
        sum = 1/sum;
        double[][] retF = scalarMultiplicationWithMatrix(FMatrix,sum);
        return retF[0];
    }
    private void takeStep(){
        double[] tT = T[head*rows*trueRow+head*trueCol+trueHead];
        int ans;
        double a = Math.random();
        double acc = 0;
        int i = 0;
        while(true){
            if(a<tT[i]+acc){
                ans = i;
                break;
            }else{
                acc+=tT[i];
                i++;
            }
        }

        trueRow = ans/(cols*head);
        trueCol = (ans%(cols*head))/4;
        trueHead = ans%4;
    }
    private double[][] scalarMultiplicationWithMatrix(double[][] A, double a){
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[0].length; j++) {
                A[i][j] *= a;
            }
        }
        return A;
    }
    private double[][] transposeMatrix(double[][] A){
        double[][] temp = new double[A[0].length][A.length];

        for(int i = 0; i < A.length; i++) {
            for(int j = 0; j < A[0].length; j++) {
                temp[j][i] = A[i][j];
            }
        }
        return temp;
    }
    private double[][] mulMatrices(double[][] A, double[][] B){
        double[][] C = new double[A.length][B[0].length];
        for(int i = 0; i < A.length; i++) {
            for(int j =  0; j < B[0].length;j++) {
                C[i][j] = 0;
                for(int m = 0; m < A[0].length; m++) {
                    C[i][j]+= A[i][m]*B[m][j];
                }
            }
        }

        return C;
    }
    private void setONothing() {
        int n = rows * cols * head;
        Onothing = new double[n];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                // System.out.println("" + i + "," + j);
                // if we are in a corner
                if (checkCornerForONothing(i, j)) {
                    // System.out.println("CORNER!");
                    setSquareToSameProbForOnothing(i, j, 0.625);
                } else {
                    int[] dis = checkDistanceToEdgeForONothing(i, j); // check distance

                    if (dis[0] == 0 && dis[1] == 1) {
                        setSquareToSameProbForOnothing(i, j, 0.5);
                    } else if (dis[0] == 0 && dis[1] >= 2) {
                        setSquareToSameProbForOnothing(i, j, 0.425);
                    } else if(dis[0] == 1 && dis[1] == 1) {
                        setSquareToSameProbForOnothing(i, j, 0.325);
                    } else if(dis[0] == 1 && dis[1] >= 2) {
                        setSquareToSameProbForOnothing(i, j, 0.25);
                    }else {
                        setSquareToSameProbForOnothing(i, j, 0.1);
                    }
                }
            }
        }

    }
    private int[] checkDistanceToEdgeForONothing(int row, int col) {
        int[] alldistances = new int[4];
        alldistances[0] = row;
        alldistances[1] = (cols - 1) - col;
        alldistances[2] = (rows - 1) - row;
        alldistances[3] = col;

        int smallestDistance = Math.min(alldistances[0], alldistances[1]);
        int nextSmallestDistance = Math.max(alldistances[0], alldistances[1]);
        for (int i = 2; i < alldistances.length; i++) {
            if(smallestDistance > Math.min(smallestDistance, alldistances[i])) {
                nextSmallestDistance = smallestDistance;
                smallestDistance = Math.min(smallestDistance, alldistances[i]);
            }else {
                nextSmallestDistance = Math.min(nextSmallestDistance, alldistances[i]);
            }

        }

        int[] distancesToReturn = {smallestDistance,nextSmallestDistance};
        return distancesToReturn;
    }
    private boolean checkCornerForONothing(int row, int col) {
        if (row == 0 && col == 0)
            return true;
        if (row == 0 && col == (cols-1))
            return true;
        if (row == (rows-1) && col == 0)
            return true;
        if (row == (rows-1) && col == (cols-1))
            return true;
        return false;
    }
    private void setSquareToSameProbForOnothing(int row, int col, double prob) {
        Onothing[col * head + row * (head * cols)] = prob;
        Onothing[col * head + row * (head * cols) + 1] = prob;
        Onothing[col * head + row * (head * cols) + 2] = prob;
        Onothing[col * head + row * (head * cols) + 3] = prob;
    }

}