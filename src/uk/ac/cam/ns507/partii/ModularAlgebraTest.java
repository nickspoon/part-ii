package uk.ac.cam.ns507.partii;

import org.jscience.mathematics.number.ModuloInteger;
import org.jscience.mathematics.number.LargeInteger;
import org.jscience.mathematics.vector.DenseMatrix;

public class ModularAlgebraTest {

    public static ModuloInteger modInt(int x) {
        return ModuloInteger.valueOf(LargeInteger.valueOf(x));
    }

    public static void main(String[] args) {
        ModuloInteger.setModulus(LargeInteger.valueOf(7));
        ModuloInteger n = modInt(3);
        System.out.println(n.plus(n));
        System.out.println(n.times(n));
        
        DenseMatrix<ModuloInteger> m = DenseMatrix.valueOf(
            new ModuloInteger[][] {
                { modInt(2), modInt(1) },
                { modInt(3), modInt(4) }
            }
        );
        System.out.println(m);
        System.out.println(m.plus(m));
        System.out.println(m.inverse());
        System.out.println(m.times(m.inverse()));
    }

}
