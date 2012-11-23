package uk.ac.cam.ns507.partii;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.jscience.mathematics.structure.Field;
import org.jscience.mathematics.vector.DenseMatrix;
import org.jscience.mathematics.vector.DenseVector;
import org.jscience.mathematics.vector.Matrix;
import org.jscience.mathematics.vector.Vector;

public class Module<K extends Field<K>> {
	private Algebra<K> lambda;
	/*
	 * structconsts[] is an array of n mxm matrices, where m is the dimension of
	 * the module and n the dimension of lambda. 
	 */
	private Matrix<K>[] structconsts;
	private int dim;
	/*
	 * ZERO is the additive identity for K.
	 */
	private K ZERO;

	public Module(Algebra<K> l, Matrix<K>[] stc, K z) {
		lambda = l;
		structconsts = stc;
		/*
		 * The dimension of the module is the number of basis elements, which is
		 * equal to the number of columns in each structure constant matrix.
		 */
		dim = stc[0].getNumberOfColumns();
		ZERO = z;
	}
	
	public Matrix<K>[] getStructureConstants() {
		return structconsts;
	}
	
	public Vector<K> computeNullSpace(Vector<K> v) {
		List<DenseVector<K>> rows = new ArrayList<DenseVector<K>>();
		for (int row = 0; row < dim; row++) {
			List<K> vect = new ArrayList<K>();
			for (int col = 0; col < lambda.getDimension(); col++) {
				K x = ZERO;
				for (int j = 0; j < dim; j++) {
					K b = v.get(j);
					x = x.plus(b.times(structconsts[col].get(j, row)));
				}
				vect.add(x);
			}
			rows.add(DenseVector.valueOf(vect));
		}
		Matrix<K> M = DenseMatrix.valueOf(rows);
		List<K> zero = Collections.nCopies(v.getDimension(), ZERO);
		DenseVector<K> z = DenseVector.valueOf(zero);
		return M.solve(z);
	}
}
