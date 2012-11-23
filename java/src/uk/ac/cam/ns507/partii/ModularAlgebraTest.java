package uk.ac.cam.ns507.partii;

import javolution.util.FastTable;
import javolution.util.Index;

import org.jscience.mathematics.number.LargeInteger;
import org.jscience.mathematics.number.ModuloInteger;
import org.jscience.mathematics.number.ModuloInteger;
import org.jscience.mathematics.vector.DenseMatrix;
import org.jscience.mathematics.vector.DenseVector;
import org.jscience.mathematics.vector.LUDecomposition;
import org.jscience.mathematics.vector.Matrix;

public class ModularAlgebraTest {

    public static ModuloInteger modInt(int x) {
        return ModuloInteger.valueOf(LargeInteger.valueOf(x));
    }
    
    public static DenseMatrix<ModuloInteger> modIntMatrix(int[][] arr) {
    	ModuloInteger[][] marr = new ModuloInteger[arr.length][arr[0].length];
    	for (int i = 0; i < arr.length; i++) {
    		for (int j = 0; j < arr[0].length; j++) {
    			marr[i][j] = modInt(arr[i][j]);
    		}
    	}
    	DenseMatrix<ModuloInteger> m = DenseMatrix.valueOf(marr);
    	return m;
    }
    
    public static void main(String[] args) {
//        ModuloInteger.setModulus(LargeInteger.valueOf(7));
//        ModuloInteger n = modInt(3);
//        System.out.println(n.plus(n));
//        System.out.println(n.times(n));
//        
//        DenseMatrix<ModuloInteger> m = DenseMatrix.valueOf(
//            new ModuloInteger[][] {
//                { modInt(2), modInt(1) },
//                { modInt(3), modInt(4) }
//            }
//        );
//        System.out.println(m);
//        System.out.println(m.plus(m));
//        System.out.println(m.inverse());
//        System.out.println(m.times(m.inverse()));
        
        ModuloInteger.setModulus(LargeInteger.valueOf(2));
        
        DenseMatrix<ModuloInteger> m = DenseMatrix.valueOf(
        		new ModuloInteger[][] {
        				{ ModuloInteger.ONE, ModuloInteger.ONE },
        				{ ModuloInteger.ONE, ModuloInteger.ONE }
        		}
        );
        DenseVector<ModuloInteger> x = DenseVector.valueOf(ModuloInteger.ONE, ModuloInteger.ONE);
        //System.out.println(m.solve(x));
        
        int[][] H1 = { {1, 0}, {0, 1} };
        int[][] H2 = { {0, 1}, {1, 0} };
        DenseMatrix<ModuloInteger> structconsts[] = new DenseMatrix[2];
        structconsts[0] = modIntMatrix(H1);
        structconsts[1] = modIntMatrix(H2); 
        Algebra<ModuloInteger> lambda = new Algebra<ModuloInteger>(structconsts);
        Module<ModuloInteger> v = new Module<ModuloInteger>(lambda, structconsts, ModuloInteger.ZERO);
        DenseVector<ModuloInteger> w = DenseVector.valueOf(ModuloInteger.ONE, ModuloInteger.ONE);
        LUDecomposition<ModuloInteger> lu = LUDecomposition.valueOf(m);
        Matrix<ModuloInteger> U = lu.getUpper(ModuloInteger.ZERO);
        Matrix<ModuloInteger> L = lu.getLower(ModuloInteger.ZERO, ModuloInteger.ONE);
        Matrix<ModuloInteger> P = lu.getPermutation(ModuloInteger.ZERO, ModuloInteger.ONE);
        System.out.println(L);
        System.out.println(U);
        System.out.println(L.times(P).times(U.times(P)));
        //System.out.println(v.computeNullSpace(w));
    }

}
