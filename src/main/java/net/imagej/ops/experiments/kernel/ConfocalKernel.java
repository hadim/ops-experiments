package net.imagej.ops.experiments.kernel;

import net.imagej.ops.OpService;
import net.imglib2.FinalDimensions;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_2D;
//import javax.media.jai.operator.ConvolveDescriptor;
//import icy.gui.dialog.MessageDialog;
//import icy.image.IcyBufferedImage;
//import icy.sequence.Sequence;
//import icy.type.DataType;

public class ConfocalKernel {
	private int _w;
	private int _h;
	private int _z;
	private double _indexImmersion;
	private double _na;
	private int _lem;
	private double _indexSpRefr;
	private double _xySampling;
	private double _zSampling;
	private double _depth;
	private double _pSize;
	private String _mName;
	private int _mObj;
	private int _lex;

	public final static String DEFAULT_MNAME = "Biorad MRC 500/600/1024";
	public final static double DEFAULT_NA = 1.4;
	public final static int DEFAULT_MOBJ = 63;
	public final static double DEFAULT_INDEXIMMERSION = 1.518;
	public final static double DEFAULT_PSIZE = 0;
	public final static int DEFAULT_LEX = 488;
	public final static int DEFAULT_LEM = 520;
	public final static double DEFAULT_INDEXSP = 1.33;
	public final static double DEFAULT_XYSAMPLING = 50.34;
	public final static double DEFAULT_ZSAMPLING = 151.37;
	public final static double DEFAULT_DEPTH = 0.0;
	public final static int DEFAULT_W = 128;
	public final static int DEFAULT_H = 128;
	public final static int DEFAULT_Z = 128;	


	ConfocalKernel() {
		setMNAME(DEFAULT_MNAME);
		setNA(DEFAULT_NA);
		setMOBJ(DEFAULT_MOBJ);
		setIndexImmersion(DEFAULT_INDEXIMMERSION);
		setPSIZE(DEFAULT_PSIZE);
		setLEX(DEFAULT_LEX);
		setLEM(DEFAULT_LEM);
		setIndexSp(DEFAULT_INDEXSP);
		setXYSAMPLING(DEFAULT_XYSAMPLING);
		setZSAMPLING(DEFAULT_ZSAMPLING);
		setDEPTH(DEFAULT_DEPTH);
		setW(DEFAULT_W);
		setH(DEFAULT_H);
		setZ(DEFAULT_Z);				
	}
	public Img compute(OpService ops){			
		double sFactor = 0.5;//Pinhole shape factor
		double mInt = 53.2;//Internal magnification		
		if(_mName == "Biorad MRC 500/600/1024")
		{
			sFactor = 0.5;
			mInt = 53.2;
			_pSize = (_pSize*sFactor)/(_mObj*mInt);
		}
		else if(_mName == "Biorad Radiance")
		{
			sFactor = 0.5;
			mInt = 73.2;
			_pSize = (_pSize*sFactor)/(_mObj*mInt);
		}
		else if(_mName == "TE2000-E C1 Head")
		{
			sFactor = 0.5;
			mInt = 1.00;
			_pSize = (_pSize*sFactor)/(_mObj*mInt);
		}
		else if(_mName == "Ti-E Perfect Focus A1R")
		{
			sFactor = 0.456;
			mInt = 1.00;		
			_pSize = (_pSize*sFactor)/(_mObj*mInt);
		}
		else if(_mName == "Olympus FV10i")
		{
			sFactor = 0.399;
			mInt = 3.8;		
			_pSize = (_pSize*sFactor)/(_mObj*mInt);
		}
		else if(_mName == "Olympus FV500")
		{
			sFactor = 0.5641896;
			mInt = 3.8;		
			_pSize = (_pSize*sFactor)/(_mObj*mInt);
		}
		else if(_mName == "Olympus FV1000")
		{
			sFactor = 0.5641896;
			mInt = 3.82;		
			_pSize = (_pSize*sFactor)/(_mObj*mInt);
		}
		else if(_mName == "Visitech Infinity")
		{
			sFactor = 0.5;
			mInt = 1.00;
			_pSize = (_pSize*sFactor)/(_mObj*mInt);
		}
		else if(_mName == "Zeiss LSM510")
		{
			sFactor = 0.5;
			mInt = 3.33;
			_pSize = (_pSize*sFactor)/(_mObj*mInt);
		}
		else if(_mName == "Zeiss LSM700")
		{
			sFactor = 0.564;
			mInt = 1.53;
			_pSize = (_pSize*sFactor)/(_mObj*mInt);
		}
		else if(_mName == "Zeiss LSM710")
		{
			sFactor = 0.564;
			mInt = 1.9048;
			_pSize = (_pSize*sFactor)/(_mObj*mInt);
		}
		else if(_mName == "Others")
		{
			System.out.println("Support for your microscope is currently unavailable. Contact us to include this.");
			return null;
		}

		final DoubleFFT_2D fft = new DoubleFFT_2D(_w, _h);
		int hc = _h/2;
		int wc = _w/2;
		int zc = _z/2;
		double kSampling = (2*Math.PI)/(_h*_xySampling); //Fourier space sampling
		double lambdaObj = _lem/_indexImmersion;//Wavelength of light inside the medium
		double lambdaSp = _lem/_indexSpRefr;//Wavelength of light inside the medium
		double k0 = (2*Math.PI)/_lem;//Wave vector
		double kObj = (2*Math.PI)/lambdaObj;//Wave vector in the immersion medium
		double kSp = (2*Math.PI)/lambdaSp;//Wave vector in the medium
		double kMax = (2*Math.PI*_na)/(_lem*kSampling);//Maximum aperture radius

		// Define the zero defocus pupil function
		Img<ComplexFloatType> pupil = ops.create().img(new FinalDimensions(_w,_h), new ComplexFloatType());// new IcyBufferedImage(_w, _h, 2, DataType.FLOAT); // channel 1 is real and channel 2 is imaginary
		float[] pupilRealBuffer = new float[_w*_h];//pupil.getDataXYAsFloat(0);//Real
		float[] pupilImagBuffer = new float[_w*_h];//pupil.getDataXYAsFloat(1);//imaginary

		//IcyBufferedImage dpupil = new IcyBufferedImage(_w, _h, 2, DataType.DOUBLE); // channel 1 is real and channel 2 is imaginary
		Img<ComplexFloatType> dpupil = ops.create().img(new FinalDimensions(_w,_h), new ComplexFloatType());
		double[] dpupilRealBuffer = new double[_w*_h];//dpupil.getDataXYAsDouble(0); //Real
		double[] dpupilImagBuffer = new double[_w*_h];//dpupil.getDataXYAsDouble(1); //imaginary

		//Calculate the cosine and the sine components
		//IcyBufferedImage ctheta = new IcyBufferedImage(_w, _h, 1, DataType.DOUBLE);
		Img<FloatType> ctheta = ops.create().img(new FinalDimensions(_w,_h), new FloatType());
		double[] cthetaBuffer = new double[_w*_h];//ctheta.getDataXYAsDouble(0);
		
		Img<FloatType> cthetaSp = ops.create().img(new FinalDimensions(_w,_h), new FloatType());
		//IcyBufferedImage cthetaSp = new IcyBufferedImage(_w, _h, 1, DataType.DOUBLE);
		double[] cthetaSpBuffer = new double[_w*_h];//cthetaSp.getDataXYAsDouble(0);
		
		Img<FloatType> stheta = ops.create().img(new FinalDimensions(_w,_h), new FloatType());
		//IcyBufferedImage stheta = new IcyBufferedImage(_w, _h, 1, DataType.DOUBLE);
		double[] sthetaBuffer = new double[_w*_h];//stheta.getDataXYAsDouble(0);
		
		Img<FloatType> sthetaSp = ops.create().img(new FinalDimensions(_w,_h), new FloatType());
		//IcyBufferedImage sthetaSp = new IcyBufferedImage(_w, _h, 1, DataType.DOUBLE);
		double[] sthetaSpBuffer = new double[_w*_h];//sthetaSp.getDataXYAsDouble(0);

		//Sequence psf3d = new Sequence();
		//psf3d.setName("Confocal PSF");
		//Sequence.phole = new Sequence();

		for (int k =  0 ; k < _z; k++)
		{// Define the defocus pupils			

			double defocus = k-zc;
			defocus = defocus*_zSampling;	

			for(int x = 0; x < _w; x++)
			{
				for(int y = 0; y < _h; y++)
				{   
					double kxy = Math.sqrt( Math.pow(x-wc, 2) + Math.pow(y-hc, 2) );

					//pupilRealBuffer[pupil.getOffset(x, y)] = ((kxy < kMax) ? 1 : 0); //Pupil bandwidth constraints
					//pupilImagBuffer[pupil.getOffset(x, y)] = 0; //Zero phase 
					pupilRealBuffer[x+y*_h] = ((kxy < kMax) ? 1 : 0); //Pupil bandwidth constraints
					pupilImagBuffer[x+y*_h] = 0; //Zero phase 

					sthetaBuffer[x + y * _h] = Math.sin( kxy * kSampling / kObj );
					sthetaBuffer[x + y * _h] = (sthetaBuffer[x + y * _h]< 0) ? 0: sthetaBuffer[x + y * _h];
					sthetaSpBuffer[x + y * _h] = Math.sin( kxy * kSampling / kSp );
					sthetaSpBuffer[x + y * _h] = (sthetaSpBuffer[x + y * _h]< 0) ? 0: sthetaSpBuffer[x + y * _h];
					cthetaBuffer[x + y * _h] = Double.MIN_VALUE + Math.sqrt(1 - Math.pow(sthetaBuffer[x + y * _h], 2));
					cthetaSpBuffer[x + y * _h] = Double.MIN_VALUE + Math.sqrt(1 - Math.pow(sthetaSpBuffer[x + y * _h], 2));
					
					dpupilRealBuffer[x + y * _h] = pupilRealBuffer[x + y * _h] * Math.cos((defocus * k0 * cthetaBuffer[x + y * _h]) + (k0 * _depth * (_indexSpRefr * cthetaSpBuffer[x + y * _h]-_indexImmersion * cthetaBuffer[x + y * _h])));
					dpupilImagBuffer[x + y * _h] = pupilRealBuffer[x + y * _h] * Math.sin((defocus * k0 * cthetaBuffer[x + y * _h]) + (k0 * _depth * (_indexSpRefr * cthetaSpBuffer[x + y * _h]-_indexImmersion * cthetaBuffer[x + y * _h])));
				}
			}
			double[] psf2d = dpupil.getDataCopyCXYAsDouble();
			fft.complexInverse(psf2d, false);

			IcyBufferedImage timg = new IcyBufferedImage(_w, _h, 1, DataType.DOUBLE);
			timg.beginUpdate();
			try{
				for(int x = 0; x < (wc+1); x++)
				{
					for(int y = 0; y < (hc+1); y++)
					{
						double val=Math.sqrt(Math.pow(psf2d[(((wc-x) + (hc-y) * _h)*2)+0],2)+Math.pow(psf2d[(((wc-x) + (hc-y) * _h)*2)+1], 2));
						timg.setDataAsDouble(x, y, 0, val);
						
						//timg.setDataAsDouble(x, y, 1, psf2d[(((wc-x) + (hc-y) * _h)*2)+1]);

					}
					for(int y = hc+1; y < _h; y++)
					{
						double val= Math.sqrt(Math.pow(psf2d[(((wc-x) + (y-hc) * _h)*2)+0], 2)+Math.pow(psf2d[(((wc-x) + (y-hc) * _h)*2)+1], 2));
						timg.setDataAsDouble(x, y, 0,val);
						//timg.setDataAsDouble(x, y, 1, psf2d[(((wc-x) + (_h+hc-y) * _h)*2)+1]);
					}

				}
				for(int x = (wc+1); x < _w; x++)
				{
					for(int y = 0; y < (hc+1); y++)
					{
						double val=Math.sqrt(Math.pow(psf2d[(((x-wc) + (hc-y) * _h)*2)+0], 2)+Math.pow(psf2d[(((x-wc) + (hc-y) * _h)*2)+1], 2));
						timg.setDataAsDouble(x, y, 0, val);
						//timg.setDataAsDouble(x, y, 1, psf2d[(((_w+wc-x) + (hc-y) * _h)*2)+1]);
					}
					for(int y = hc+1; y < _h; y++)
					{
						double val=Math.sqrt(Math.pow(psf2d[(((x-wc) + (y-hc) * _h)*2)+0], 2)+Math.pow(psf2d[(((x-wc) + (y-hc) * _h)*2)+1],2));
						timg.setDataAsDouble(x, y, 0, val);
						//timg.setDataAsDouble(x, y, 1, psf2d[(((_w+wc-x) + (_h+hc-y) * _h)*2)+1]);
					}
				}

			}finally {
				timg.endUpdate();
			}

			psf3d.addImage(timg);
		}

		return psf3d;

	}
	public int getW() {
		return _w;
	}
	public int getH() {
		return _h;
	}
	public int getZ() {
		return _z;
	}
	public double getIndexRefr() {
		return _indexImmersion;
	}
	public double getNA() {
		return _na;
	}
	public int getLEM() {
		return _lem;
	}
	public double getSA() {
		return _indexSpRefr;
	}
	public double getXYSAMPLING() {
		return _xySampling;
	}
	public double getZSAMPLING() {
		return _zSampling;
	}
	public double getDEPTH() {
		return _depth;
	}

	public void setW(int src) {
		_w = src;
	}
	public void setH(int src) {
		_h = src;
	}
	public void setZ(int src) {
		_z = src;
	}
	public void setIndexImmersion(double src) {
		_indexImmersion = src;
	}
	public void setNA(double src) {
		_na = src;
	}
	public void setLEM(int src) {
		_lem = src;
	}
	public void setIndexSp(double src) {
		_indexSpRefr = src;
	}
	public void setXYSAMPLING(double src) {
		_xySampling = src;
	}
	public void setZSAMPLING(double src) {
		_zSampling = src;
	}
	public void setDEPTH(double src) {
		_depth = src;
	}
	public void setPSIZE(double src) {
		_pSize = src;
	}
	public void setMNAME(String src) {
		_mName = src;
	}
	public void setMOBJ(int src) {
		_mObj = src;
	}
	public void setLEX(int src) {
		_lex = src;
	}	
}
