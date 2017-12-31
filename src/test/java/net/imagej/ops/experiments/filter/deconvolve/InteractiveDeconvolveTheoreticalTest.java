package net.imagej.ops.experiments.filter.deconvolve;

import java.io.IOException;

import net.imagej.ImageJ;
import net.imagej.ops.filter.pad.DefaultPadInputFFT;
import net.imglib2.FinalDimensions;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.outofbounds.OutOfBoundsMirrorFactory;
import net.imglib2.outofbounds.OutOfBoundsMirrorFactory.Boundary;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;

public class InteractiveDeconvolveTheoreticalTest<T extends RealType<T> & NativeType<T>> {

	final static ImageJ ij = new ImageJ();

	public static <T extends RealType<T> & NativeType<T>> void main(final String[] args) throws IOException {

		String libPathProperty = System.getProperty("java.library.path");
		System.out.println("Lib path:" + libPathProperty);

		ij.launch(args);

		String inputName = "C:/Users/bnorthan/Dropbox/Deconvolution_Test_Set/FromJeanYves/Slide_17015-02_512_2.tif";
		String psfName = "../ops-images/deconvolution/PSF-CElegans-CY3-cropped.tif";

		@SuppressWarnings("unchecked")
		Img<T> img = (Img<T>) ij.dataset().open(inputName).getImgPlus().getImg();
		Img<FloatType> imgF = ij.op().convert().float32(img);

		// create the diffraction based image 
		Img<FloatType> psf = (Img) ij.op().create().kernelDiffraction(new FinalDimensions(64, 64, 50), 1.4, 550E-09,
				1.3, 1.5, 162.5E-9, 280E-9, 0, new FloatType());

		// normalize PSF
		/*FloatType sum = new FloatType(ij.op().stats().sum(psf).getRealFloat());
		psf = (Img<FloatType>) ij.op().math().divide(psf, sum);

		// extend image
		RandomAccessibleInterval<FloatType> extendedImage = (RandomAccessibleInterval<FloatType>) ij.op().run(
				DefaultPadInputFFT.class, imgF,
				new FinalDimensions(img.dimension(0), img.dimension(1), img.dimension(2)), true,
				new OutOfBoundsMirrorFactory<>(Boundary.SINGLE));*/

		ij.ui().show("bars ", img);
		ij.ui().show("psf", psf);

		int iterations = 10;
		int pad = 0;

		long startTime, endTime;

		// run Ops Richardson Lucy

		startTime = System.currentTimeMillis();

		Img<FloatType> deconvolved = (Img<FloatType>) ij.op().deconvolve().richardsonLucyTV(imgF, psf, null, null, null,
				null, null, iterations, true, true, 0.0002f);
		
		ij.ui().show("deconvolved", deconvolved);
		/*
		 * Img<FloatType> deconvolved = (Img<FloatType>)
		 * ij.op().deconvolve().richardsonLucy(imgF, psf, new long[] { pad, pad,
		 * pad }, null, null, null, null ,10, true, true);
		 * 
		 * 
		 * endTime = System.currentTimeMillis();
		 * 
		 * ij.log().info("Total execution time (Ops) is: " + (endTime -
		 * startTime));
		 * 
		 * ij.ui().show("Richardson Lucy deconvolved", deconvolved);
		 * 
		 * // run Cuda Richardson Lucy op
		 * 
		 * startTime = System.currentTimeMillis();
		 * 
		 * RandomAccessibleInterval<FloatType> outputCuda =
		 * (RandomAccessibleInterval<FloatType>) ij.op()
		 * .run(CudaRichardsonLucyOp.class, imgF, psf, new long[] { pad, pad,
		 * pad }, iterations);
		 * 
		 * endTime = System.currentTimeMillis();
		 * 
		 * ij.log().info("Total execution time (Cuda) is: " + (endTime -
		 * startTime));
		 * 
		 * ij.ui().show("cuda op deconvolved", outputCuda);
		 * 
		 * // run MKL Richardson Lucy
		 * 
		 * startTime = System.currentTimeMillis();
		 * 
		 * // run MKL Richardson Lucy op RandomAccessibleInterval<FloatType>
		 * outputMKL = (RandomAccessibleInterval<FloatType>) ij.op()
		 * .run(MKLRichardsonLucyOp.class, imgF, psf, new long[] { pad, pad, pad
		 * }, iterations);
		 * 
		 * endTime = System.currentTimeMillis();
		 * 
		 * ij.log().info("Total execution time (MKL) is: " + (endTime -
		 * startTime));
		 * 
		 * ij.ui().show("mkl op deconvolved", outputMKL);
		 */

	}

}
