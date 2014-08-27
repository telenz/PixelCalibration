PixelCalibration
================

source code needed for calculating calibration factors

	mkdir PixelCalibration
	cd PixelCalibration
	cmsrel CMSSW_5_3_15
	cd CMSSW_5_3_15
	git clone https://github.com/telenz/PixelCalibration.git src
	cd src
	#git checkout tags/T_PIXCALIB_0
	cmsenv
	git cms-cvs-history import CMSSW_5_3_11 FWCore/Version
	scram b -j24
	#end