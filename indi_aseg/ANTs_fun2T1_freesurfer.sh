# brainmask.nii  label- label_area.nii
betfundir=/mnt/d/jingzongfmri/data/Indi/fmri_indi_v1_new/Betindicoreg
braindir=/mnt/d/jingzongfmri/data/Indi/fmri_indi_v1_new/T1Indi
fundir=/mnt/d/jingzongfmri/data/Indi/fmri_indi_v1_new/FunImgAR

cd $betfundir
for file in $(ls)
	do
		cd $fundir/$file
		fpath=($betfundir/$file)
		filef=$(ls $fpath/Bet*)
		betfun=$filef

		tpath=($braindir/$file)
		filet=$(ls $tpath/RS_brain*)
		fileT=$filet

		
		#antsRegistrationSyN.sh -d 3 -f $fileT -m $betfun -o reg -t 'r'
		ffpath=($fundir/$file)
		fileff=$(ls $ffpath/*4D*)
		funff=$fileff
		#antsApplyTransforms -i $funff -o $file"_rsfMRI_reg.nii" -r $fileT -t reg0GenericAffine.mat -e 3

                                
	done


