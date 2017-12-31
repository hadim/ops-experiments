package net.imagej.experiments.transform;

import java.io.IOException;

import net.imagej.ImageJ;
import net.imagej.ops.special.computer.BinaryComputerOp;
import net.imagej.ops.special.computer.Computers;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.interpolation.randomaccess.LanczosInterpolatorFactory;
import net.imglib2.realtransform.AffineTransform;
import net.imglib2.realtransform.AffineTransform2D;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

public class InteractiveTransformOpTest<T extends RealType<T> & NativeType<T>> {

	static String workingDir = "./images/";

	static String imName = workingDir + "lowresbridge.tif";

	public static <T extends RealType<T> & NativeType<T>> void main(final String[] args) throws IOException {

		ImageJ ij;

		ij = net.imagej.Main.launch();

		Img<T> image = (Img<T>) (ij.dataset().open(imName).getImgPlus().getImg());

		Img<T> imageTransformed = ij.op().create().img(image);

		ij.ui().show("Original", image);

		BinaryComputerOp<RandomAccessibleInterval<T>, AffineTransform2D, IterableInterval<T>> transformOp = null;

		// define a affine transform
		AffineTransform2D trans = new AffineTransform2D();
/*
		trans.translate(-5., 0);

		transformOp.compute(image, trans, imageTransformed);

		ij.ui().show("shifted", imageTransformed);

		trans = new AffineTransform2D();
		trans.rotate(1);

		transformOp.compute(image, trans, imageTransformed);

		ij.ui().show("rotated origin", imageTransformed);
		
		trans = new AffineTransform2D();
		
		trans.translate(-image.dimension(0)/2,-image.dimension(0)/2);
		trans.scale(1.5);
		trans.rotate(1);
		trans.translate(image.dimension(0)/2,image.dimension(0)/2);
		
		transformOp.compute(image, trans, imageTransformed);
		
		ij.ui().show("rotated center 1", imageTransformed);
		
		
		//ij.op().transform().realT
		
*/
		//trans.translate(-image.dimension(0)/2,-image.dimension(0)/2);
		//trans.scale(1.5);
		trans.translate(-image.dimension(0)/2,-image.dimension(0)/2);
		trans.rotate(1);
		trans.scale(0.5);
		trans.translate(image.dimension(0)/2,image.dimension(0)/2);
		
		

		ij.ui().show("rotated center 2", ij.op().transform().realTransform(image, trans));
		
		// define a affine transform
		trans = new AffineTransform2D();
		
		double[][] mat=new double[][]{new double[]{1,0.45,0}, new double[]{0,1.,0},new double[]{0,0,1}};
		trans.set(mat);
		
		Interval shearInterval=new FinalInterval(new long[]{0,0},new long[]{2*image.dimension(0)-1,image.dimension(1)-1});
		
		ij.ui().show("shear", ij.op().transform().realTransform(image, trans, shearInterval));
		
		
	/*	
		trans = new AffineTransform2D();
		trans.scale(.5);
		
		transformOp.compute(image, trans, imageTransformed);
		
		ij.ui().show("scaled", imageTransformed);*/

	}

}
