PixelCalibration
================

source code needed for calculating calibration factors

	mkdir PixelCalibration
	cd PixelCalibration
	cmsrel CMSSW_5_3_15
	cd CMSSW_5_3_15
	git clone https://github.com/telenz/PixelCalibration.git src
	#git checkout tags/T_PIXCALIB_0
	cd src
	cmsenv
	git cms-cvs-history import CMSSW_5_3_11 FWCore/Version
	scram b -j24
	#end