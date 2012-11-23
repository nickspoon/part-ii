package uk.ac.cam.ns507.partii;

import org.jscience.mathematics.structure.Field;
import org.jscience.mathematics.vector.Matrix;

public class Algebra<K extends Field<K>> {
	/*
	 * structconsts[] is an array of n nxn matrices, where n is the dimension of
	 * the algebra.
	 */
	private Matrix<K> structconsts[];
	private int dim;

	public Algebra(Matrix<K>[] stc) {
		structconsts = stc;
		dim = stc[0].getNumberOfColumns();
	}

	public Matrix<K>[] getStructureConstants() {
		return structconsts;
	}
	
	public int getDimension() {
		return dim;
	}
}
