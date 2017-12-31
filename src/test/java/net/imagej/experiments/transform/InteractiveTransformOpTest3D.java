package net.imagej.experiments.transform;

import java.io.IOException;

import net.imagej.ImageJ;
import net.imagej.ops.special.computer.BinaryComputerOp;
import net.imagej.ops.special.computer.Computers;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.interpolation.InterpolatorFactory;
import net.imglib2.interpolation.randomaccess.LanczosInterpolatorFactory;
import net.imglib2.outofbounds.OutOfBoundsConstantValueFactory;
import net.imglib2.outofbounds.OutOfBoundsFactory;
import net.imglib2.realtransform.AffineGet;
import net.imglib2.realtransform.AffineTransform;
import net.imglib2.realtransform.AffineTransform2D;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.RealViews;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.util.Util;
import net.imglib2.view.Views;

public class InteractiveTransformOpTest3D<T extends RealType<T> & NativeType<T>> {

	static String workingDir = "../ops-images/deconvolution/";
	
	static String imName = workingDir + "Bars-G10-P15-stack.tif";
	
	public static <T extends RealType<T> & NativeType<T>> void main(final String[] args) throws IOException {

		ImageJ ij;

		ij = net.imagej.Main.launch();

		Img<T> image = (Img<T>) (ij.dataset().open(imName).getImgPlus().getImg());

		Img<T> imageTransformed = ij.op().create().img(image);

		ij.ui().show("Original", image);
		
		AffineTransform3D transform=new AffineTransform3D();
		
		T zero = Util.getTypeFromInterval(image).copy();
		zero.setZero();
		OutOfBoundsFactory<T, RandomAccessibleInterval<T>> outOfBoundsFactory = new OutOfBoundsConstantValueFactory<>(zero);
		
		InterpolatorFactory<T, RandomAccessible<T>> interpolator =
				new LanczosInterpolatorFactory<>();
		
		double[][] mat=new double[][]{new double[]{1,0.45,0,0}, new double[]{0,1.,0,0},new double[]{0,0,1,0},new double[]{0,0,0,1}};
		transform.set(mat);
		
		Interval outputInterval=new FinalInterval(new long[]{0,0,0},new long[]{2*image.dimension(0)-1,image.dimension(1)-1,image.dimension(1)-1});
		
		Views.interval(Views.raster(RealViews.affineReal(Views.interpolate(
				Views.extend(image, outOfBoundsFactory), interpolator),
				(AffineGet) transform)), outputInterval);
		
		//ij.ui().show("shear", ij.op().transform().realTransform(image, transform));//, shearInterval));



	}

}
