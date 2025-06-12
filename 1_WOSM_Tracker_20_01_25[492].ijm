/* Simple Object-Tracking Program

	WOSM_tracker.txt

	Justin Molloy
	Centre for Mechanochemical Cell Biology (L1.01)
	Warwick Medical School
	Coventry CV4 7AL
	M: +44 (0) 7986143550
	E: justin.molloy@warwick.ac.uk
	W: mechanochemistry.org

	18-Nov-2022
	Modified: 06-Nov-2024
		  20-Jan-2025 -- output immediate velocity with x, y coordinates

ver3:
	The program identifies spots using ImageJ "FindMaxima" function
	You can find bright or dark spots (or both).
	The program then tracks spot locations between frames within a limiting radius.
	Trims out short tracks.
	Reorders and renumbers remaining tracks
	Performs fitting on MSD vs. dT plots 
	Plots all tracks onto the stack overlay for visual inspection
	Saves all data to the "Results" window

ToDo: 
	Discard immobile spots
	Compute super-resolution fitting for each tracked spot
	Compute persistence of motion
*/

requires("1.54b");
call("java.lang.System.gc");

// create some dinky colours:
// see:   https://en.wikipedia.org/wiki/List_of_colors
avocado= "#568203"; lime =   "#D0FF14"; azure =  "#007FFF"; beau =   "#BCD4E6";
capri=   "#00BFFF"; cyan =   "#00B7EB"; orchid = "#E29CD2"; paleOrch ="#FAE6FA";
brown =  "#88540B"; copper = "#DA8A67"; imperial="#ED2939"; pink =   "#EFBBCC";
colors = newArray(avocado, lime, azure, beau, capri, cyan, orchid, paleOrch, brown, copper);
colors = Array.concat(colors, imperial, pink,"blue", "green", "magenta", "orange", "yellow","red","pink");

var um = fromCharCode(181)+"m";	var um2= um+fromCharCode(178);
var umPerSec= um+"/s";	var um2PerSec=um2+"/s";

// create a Laplacian of a Gaussian kernal for spot-finding
var LoG7 = "text1=[";
LoG7 +=" 0  0 -1 -1 -1  0  0\n";
LoG7 +=" 0 -1 -1 -1 -1 -1  0\n";
LoG7 +="-1 -1  2  5  2 -1 -1\n";
LoG7 +="-1 -1  5 10  5 -1 -1\n";
LoG7 +="-1 -1  2  5  2 -1 -1\n";
LoG7 +=" 0 -1 -1 -1 -1 -1  0\n";
LoG7 +=" 0  0 -1 -1 -1  0  0 ]";


//Initialise things
//===========
run("Set Measurements...", "display redirect=None decimal=3");
if (isOpen("Plot of Results")){
	selectWindow("Plot of Results");
	close();
}
if (isOpen("MSD vs dT STACK")){
	selectWindow("MSD vs dT STACK");
	close();
}

// save MSDdt and XY coordinates to ".CSV" files
//==============================================
saveXY = false;
saveMSD = false;
saveResults =false;
path = getDir("image");
saveXYfile = path+"XY_coord.csv";
saveBrownian_MSD_file = path+"Brownian_MSD_dT.csv";
saveMotor_MSD_file = path+"Motor_MSD_dT.csv";
saveResultsfile = path+"Results.csv";

// Load the Ome.tif Metadata
print(LoadMetaData());

// plotting limits for MSD vs dT
//=====================
var plotdTmax=2; var plotXmax=2; var plotYmax=2;

// Data Fitting rage and minimum acceptable R-value for fits

var alpha = 1.6; // exponential discrimator for Brownian vs motorised
var minR = 0.5;
var fitRange= 80;
 // set to 80%

// smoothing window size for Lp analysis
var sm = 9;

//User dialog 
//===========
// defaults...
var mintracklen = 10;
var radius = 10;
var dark = true;
var subPrev = false;
var bullseye = false;
var trackID = true;
var filter = false;
var FFTfilter = false;
var nObs = 100;
var lineThick = 2;
var trackFontSize=22;


// Put up the initial Dialog box
Dialog.create("Auto Detect");
Dialog.addNumber("How many objects:", nObs);
Dialog.addCheckbox("Filter (FFT_Bandpass)",FFTfilter);
Dialog.addCheckbox("Subtract Previous Frame",subPrev);
Dialog.addCheckbox("Filter (LoG)+GaussBlur",filter);
Dialog.addCheckbox("Dark Spots",dark);

Dialog.show();

nObs=Dialog.getNumber();
FFTfilter= Dialog.getCheckbox();
subPrev = Dialog.getCheckbox();
filter= Dialog.getCheckbox();
dark = Dialog.getCheckbox();

// Set the FindMaxima Option.. "dark spots" require "light" setting!
opt= " ";
if (dark) opt= " light "; 

// Filter the stack as specified
if (filter || FFTfilter) run("Duplicate...", "duplicate");
if (FFTfilter) run("Bandpass Filter...", "filter_large=40 filter_small=3 suppress=None tolerance=5 process");

// subPrevFrame

if (subPrev) {
	imageTitle=getTitle();
	run("Duplicate...", "title=shifted duplicate range=2-"+ d2s(nSlices(),0) );
	selectWindow(imageTitle);
	setSlice(nSlices());
	run("Delete Slice");
	imageCalculator("Subtract create 32-bit stack","shifted",imageTitle);
	rename("temp");
	selectWindow(imageTitle);
	close();
	selectWindow("shifted");
	close();
	selectWindow("temp");
	rename(imageTitle);
}

if (filter) {
	run("Convolve...", LoG7+" normalize stack");
	run("Gaussian Blur...", "sigma=2 stack");
}


// get info about selected stack
//======================
inputStack = getTitle();
//run("Remove Overlay");
nFrames = nSlices();
var xCal; var yCal; var zCal;
getVoxelSize(xCal, yCal, zCal, unit);
if (Stack.getFrameInterval() >0) zCal = Stack.getFrameInterval();
//print ("Reset Frame interval to ",zCal);
setVoxelSize(xCal, yCal, zCal, unit);

// find the correct prominence value based on expected number of spots
Prominence = autoProminence(nObs,opt);

// Now display the Main User Dialog Box
Dialog.create("Find Maxima");
Dialog.setInsets(0, 50, 0); Dialog.addMessage("Tracking parameters:",16,"red");
Dialog.addNumber("Suggested Prominence:",Prominence);
Dialog.addNumber("Min Track length:", mintracklen);
Dialog.addNumber("Max link search radius:", radius);

Dialog.setInsets(0, 50, 0); Dialog.addMessage("MSD vs dT Fitting:",16,"red");
Dialog.addNumber("Min exponent (alpha) for motors", alpha);
Dialog.addNumber("Minimum R-value for plots", minR);
Dialog.addNumber("Limit range for fitting 10-100 (%)", fitRange);

Dialog.setInsets(0, 50, 0);Dialog.addMessage("Output settings:",16,"red");
Dialog.addNumber("Track Line Thickness", lineThick);
Dialog.addNumber("Track Font Size", trackFontSize);
Dialog.setInsets(10, 20, 0);Dialog.addCheckbox("Plot Bullseyes", bullseye);
Dialog.addCheckbox("Plot Track IDs", trackID);
Dialog.addCheckbox("Save \"ID, x, y\" as Excel file", saveXY);
Dialog.addCheckbox("Save \"ID, dT, MSD, err\" as Excel file", saveMSD);
Dialog.addCheckbox("Save \"Results Table\" as Excel file", saveResults);

Dialog.show();

//Get User input values
Prominence = Dialog.getNumber();
mintracklen = Dialog.getNumber();
radius = Dialog.getNumber();

alpha = Dialog.getNumber();
minR = Dialog.getNumber();
fitRange = Dialog.getNumber();

lineThick = Dialog.getNumber();
trackFontSize = Dialog.getNumber();
bullseye = Dialog.getCheckbox();
trackID = Dialog.getCheckbox();
saveXY = Dialog.getCheckbox();
saveMSD = Dialog.getCheckbox();
saveResults = Dialog.getCheckbox();

// check the values are "in range" and valid
if (fitRange<10) fitRange=10;
if (fitRange>100) fitRange=100;
if (minR<0) minR=0;
if (minR>1) minR=1;
if ((saveXY) && File.exists(saveXYfile)) File.delete(saveXYfile);
if ((saveMSD) && File.exists(saveBrownian_MSD_file)) File.delete(saveBrownian_MSD_file);
if ((saveMSD) && File.exists(saveMotor_MSD_file)) File.delete(saveMotor_MSD_file);

if ((saveResults) && File.exists(saveResultsfile)) File.delete(saveResultsfile);

//================== Find Spots and track them ==================
setBatchMode("hide");

x_vals=newArray();
y_vals=newArray();
slice_num=newArray();
track_num=newArray();
track_num_pts=newArray();
trackpts=newArray();

showStatus("Finding Spots and Linking..");
for (i=1; i<=nFrames; i++) {		//<-- Outer loop
	showProgress(i, nFrames);
	selectWindow(inputStack);
	setSlice(i);

	ct=0;
	xpts=newArray();
	ypts=newArray();
	nslice=newArray();
	nTrack=newArray();

// Find spots
//===========
	if (Prominence > 0){
		run("Find Maxima...", "prominence="+ Prominence + opt +" output=List");
		dk=nResults;
		for (j=0; j<dk; j+=1){
			xpts[ct]=getResult("X",j);
			ypts[ct]=getResult("Y",j);
			nslice[ct]=i;
			nTrack[ct]=ct;
			ct+=1;
		}
	}

// Note: This algorithm needs to be modifed.. it prioritises 
// matching new spots to old. It should prioritise matching 
// old spots with new ones!
//
// Link spots to nearest neighbour on the previous frame
//======================================================
// First sort the arrays to speed-up linking search.

	Array.sort(xpts, ypts);
	
// find the neighbour spot on current frame compared to previous frame
// This code needs to be re-written to prioritise finding olds spots to match the new spots 
// not as it is being done currently.. the reason this matters is we want the search radius
// originate from the new spot positions.. not the old ones.

	if(i>1){
		for (j=0; j<ct ;j+=1){
			distmin=radius;
			k=0; xmax=0; killIndex=-1;
			while ((k<oldct) && (xmax<=radius)){
				// oldct = number of tracks on previous frame
				distx=oldxpts[k]-xpts[j];
						// xmax is the last x value tested. Since xpts are ranked 
				disty=oldypts[k]-ypts[j];
						// stop searching the list when xmax > radius
				dist=sqrt(distx*distx + disty*disty);
			// Euclidian distance between current (jth) track point and 
				if (dist<distmin){
								// the kth point found on in the previous frame
					distmin=dist;
								// Euclidian minimum distance
					nTrack[j] = oldnTrack[k];
					killIndex=k;
								// the array index of the nearest point
				}
				xmax=maxOf(distx,xmax);
							// keep searching until the next x value 
				k+=1;
											// exceeds the search radius limit 
			}
	
			if (killIndex == -1){
								// no near-neighbour object found
				nTrack[j]= nextTrack;
				trackpts[nTrack[j]]=1;
				nextTrack+=1;
			} else {
				trackpts[nTrack[j]]+=1;
							// we found a linked point remove it from search list
				oldnTrack=Array.deleteIndex(oldnTrack, killIndex);
				oldxpts=Array.deleteIndex(oldxpts, killIndex);
				oldypts=Array.deleteIndex(oldypts, killIndex);
				oldct=oldnTrack.length;
							// update the "while loop" counter
			}
		}

	} else {
													// Special case.. first frame
		nextTrack=ct;
		for (j=0;j<nextTrack;j+=1) trackpts[nTrack[j]]=1;
	}

// "new spots" are now "old spots" 
//================================
	oldxpts=Array.copy(xpts);
	oldypts=Array.copy(ypts);
	oldnTrack=Array.copy(nTrack);
	oldct=nTrack.length;

// concatenate the data into mega arrays
//======================================
	slice_num= Array.concat(slice_num, nslice);
	track_num= Array.concat(track_num, nTrack);
	x_vals= Array.concat(x_vals,xpts);
	y_vals= Array.concat(y_vals,ypts);
	

showProgress(i,nFrames);

} // <--- outer loop.... do next slice!!!

//=========== finished spot detection and linking ===========

// Remove short tracks
//====================

showStatus("Removing short tracks");

for (i=0; i<slice_num.length; i+=1){
	track_num_pts[i]=trackpts[track_num[i]];
}

Array.sort(track_num_pts, track_num, slice_num, x_vals, y_vals);

ct=0;
while ((ct<slice_num.length) && (track_num_pts[ct] < mintracklen)) ct+=1;

slice_num= Array.slice(slice_num, ct, slice_num.length);
track_num= Array.slice(track_num, ct, track_num.length);
x_vals= Array.slice(x_vals, ct, x_vals.length);
y_vals= Array.slice(y_vals, ct, y_vals.length);
track_num_pts=Array.slice(track_num_pts, ct, track_num_pts.length);

// Re-Number
//==========
Array.sort(track_num, track_num_pts, slice_num, x_vals, y_vals);

// we can now reuse the "trackpts" array
trackpts= newArray();
trackIndex= newArray();
ct=0; trackIndex[ct]=0;

for (i=0;i<track_num.length-1; i+=1){
	oldnum=track_num[i];
	track_num[i]=ct;
	trackpts[ct]=track_num_pts[i];
	trackIndex[ct+1]= trackIndex[ct] + trackpts[ct];
	if (track_num[i+1]!= oldnum) ct+=1;
}

track_num[track_num.length-1]=ct;
nTracks=trackpts.length;


// Compute smoothed trajectories
//----------------------------------------------------
// Smooth the xval and yval data over a "sm"-point running window.
// sm = window size
	sm2 = sm/2;
	smoothx = Array.copy(x_vals);
	smoothy = Array.copy(y_vals);
	st=0;

	for (i=0; i<nTracks; i+=1){			// iterate through all tracks
		nd = trackpts[i]+st-sm;			// st and nd points for each track
		for (j=st; j<nd; j+=1){
			accumx=0;
			accumy=0;
			for (k=0; k<sm; k++){
				accumx+= x_vals[k+j];		// sum over sm points
				accumy+= y_vals[k+j];
			}
			smoothx[j+sm2] = accumx/sm;		// find average
			smoothy[j+sm2] = accumy/sm;	
		}		
	st+= trackpts[i];						// move to new start point in arrays
	}

// optionally replace the tracks with smoothed values
//x_vals= Array.copy(smoothx);
//y_vals= Array.copy(smoothy);

// Compute MSD vs dT
//==================
st=0; // starting offset into x_vals and y_vals arrays
dT=newArray(slice_num.length);
MSD=newArray(slice_num.length);
MSD=Array.fill(MSD,0);
ndTs=Array.copy(MSD);
err=Array.copy(MSD);
weights=Array.copy(MSD);

showStatus("Computing MSD vs dTs");
// iterate across all timepoint pairs "(k) and (k+d)" and all tracks (j) 
// Three nested loops so this could take a while.... In fact, it's surprisingly quick!
for (i=0; i<nTracks; i+=1){
	nPts=trackpts[i];
	for (j=1; j<nPts; j+=1){
		for (k=st; k<(st+nPts-j); k+=1){
			MSD[st+j] = MSD[st+j] + ( (x_vals[k+j]-x_vals[k])*(x_vals[k+j]-x_vals[k]) ); 
			MSD[st+j] = MSD[st+j] + ( (y_vals[k+j]-y_vals[k])*(y_vals[k+j]-y_vals[k]) );
			ndTs[st+j]+=1;
		}
		dT[st+j]=j*zCal;
	}
	st+=nPts;
showProgress(i,nTracks);
}

// Compute the errors and Fit-weightings for the MSD estimates (based on number of samples)
for (i=0; i<MSD.length; i+=1){
	if (ndTs[i]>0){
		MSD[i]=(MSD[i]/ndTs[i])*xCal*yCal;	
		err[i] = MSD[i]/ndTs[i];
		weights[i] = 1 /(err[i] * err[i]);
	}
}

/* First-Pass compute the following:
   ==================================
   where 	i = track counter 
   and 		j = point counter within track
   meanVel = Average "immediate" velocity
   ==================================
*/
        immedVels = newArray(x_vals.length);
	meanVel = newArray(nTracks);
	Array.fill(meanVel, 0);
	sigmadsk = Array.copy(meanVel);

st=0;
	for (i=0; i<nTracks; i+=1){			// iterate through all tracks
		ct=0;
		nd = trackpts[i]+st-1;			// st and nd points for each track
		Theta1=0;
		for (j=st; j<nd; j+=1){
			dx=x_vals[j+1]-x_vals[j];
			dy=y_vals[j+1]-y_vals[j];
			chord2 = dx*dx + dy*dy;
			dist=0;
			if (chord2>0) dist=sqrt(chord2)*xCal;
			sigmadsk[i]+= dist;
			immedVels[j]= dist/zCal;
			ct+=1;
		} // do next point

	meanVel[i] = sigmadsk[i] /(zCal*ct);		// units of um/second
	st+= trackpts[i];
	} // do next track
	
//================================================
// Output selected data to .csv (EXCEL) data files
//================================================

if (saveXY){
showStatus("Output x,y as csv file");
	st=0;
	for (i=0; i<nTracks; i+=1){
		nPts=trackpts[i];
		for (j=st; j<(st+nPts); j+=1){
			str= d2s(i,0) +","+ d2s(x_vals[j]*xCal,3) +","+ d2s(y_vals[j]*yCal,3)+","+ d2s(immedVels[j],3);
			File.append(str, saveXYfile);
		}
		File.append(" ", saveXYfile);
		st+=nPts;
		showProgress(i, nTracks);
	}
}


setBatchMode(false);

//plot tracks on movie
//====================
	plotOverlays();

//create MSD vs dT plot movie
//===========================
	plotMSD_dT_Movie(0,nTracks);

//Save Results.csv to disc
//========================
	saveAs("Results", saveResultsfile);

// collect memory garbage
//=======================
	call("java.lang.System.gc");

//====================================== FINISHED ================================



// Plotting functions below here....................
//==========================================
function plotMSD_dT_Movie(trackst, tracknd) {
//==========================================
xLab= "dT(s)"; yLab= "MSD (" + um2 + ")";
run("Clear Results");
brownWalk=false;
motor=false;

sliceNo=1;
motorStr="";
nVel=0;
nDlat=0;
nS=0;

meanBrownianMSD=newArray(25);
Array.fill(meanBrownianMSD,0);
meanMotorMSD = Array.copy(meanBrownianMSD);
motorCt= Array.copy(meanBrownianMSD);
brownCt= Array.copy(meanBrownianMSD);

xpts= newArray();
ypts= newArray();
Vel= newArray();
Dlat= newArray();
S_exp= newArray();

showStatus("Plotting MSD vs dT graphs");

Plot.create("Dynamic Plots", xLab, yLab);
Plot.setFrameSize(512,512);
//Plot.setLimits(0, plotXmax, 0, plotYmax);
Plot.setFontSize(16, "normal");
Plot.setLineWidth(2);

ct=0;
for (j=trackst; j<tracknd; j+=1){
	st= trackIndex[j];
	nd= trackIndex[j]+trackpts[j];

//	range= minOf(trackpts[j]/2, plotdTmax/zCal);
	range= trackpts[j]*(fitRange/100);
	xpts= newArray(); xpts=Array.slice(dT,st,nd);
	ypts= newArray(); ypts=Array.slice(MSD,st,nd);
	errbar= newArray(); errbar=Array.slice(err,st,nd);
	wts= newArray(); wts=Array.slice(weights,st,nd);
	Array.getStatistics(ypts, minY, maxY, mean, stdDev);

// Fit to a parabola (motorized movement)
	Fit.doWeightedFit("y = a*(x*x) + b", Array.slice(xpts,0,range) , Array.slice(ypts,0,range), Array.slice(wts,0,range));
	Vfit = newArray();
	for (i=0; i<range; i+=1) Vfit[i] = Fit.f(xpts[i]);
	Vel[ct] = sqrt(Fit.p(0)); 
	VRsq = Fit.rSquared;
	
// Fit to a Brownian Walk
	Fit.doWeightedFit("y = a*x + b", Array.slice(xpts,0,range) , Array.slice(ypts,0,range), Array.slice(wts,0,range));
	Dfit = newArray();
	for (i=0; i<range; i+=1) Dfit[i] = Fit.f(xpts[i]);

// 2-dimensional Dlat = gradient/4
	Dlat[ct] = Fit.p(0)/4; 
	DlatRsq = Fit.rSquared;

// Best fit to an arbitrary power-law (with offset)
	Fit.doWeightedFit("y = a*pow(x,b) + c", Array.slice(xpts,0,range) , Array.slice(ypts,0,range), Array.slice(wts,0,range));
	Sfit = newArray();
	for (i=0; i<range; i+=1) Sfit[i] = Fit.f(xpts[i]);

	S_exp[ct] = Fit.p(1); 
	S_expRsq = Fit.rSquared;
	
// minR is the minimum accepatable fit value

valid=false;
//	if ((VRsq>DlatRsq) && (VRsq>minR)) {
	if ((S_exp[ct] >= alpha) && (S_expRsq>minR)) {
		setResult("ID_Vel", nVel, j);
		setResult("Vel", nVel, Vel[ct]);
		setResult("R_Vel", nVel, VRsq);
		setResult("n_Vel", nVel, trackpts[j]);
		nVel+=1;
		valid=true;
	}
//	if ((DlatRsq>VRsq) && (DlatRsq>minR)) { 
	if ((S_exp[ct] < alpha) && (S_expRsq>minR)) {
		setResult("ID_Dlat", nDlat, j);
		setResult("Dlat", nDlat, Dlat[ct]);
		setResult("R_Dlat", nDlat, DlatRsq);
		setResult("n_Dlat", nDlat,trackpts[j]);
		nDlat+=1;
		valid=true;
	}
	
	if (S_expRsq>minR) { 
		setResult("ID_Exp", nS, j);
		setResult("Exp", nS, S_exp[ct]);
		setResult("R_Exp", nS, S_expRsq);
		setResult("Immed_Vel", nS, meanVel[j]);
		setResult("n_Exp", nS, trackpts[j]);
		nS+=1;
		valid=true;
	}

//Plot the graphs and create a list of motorized spots 
 if (valid) {
	Plot.setColor(colors[j % colors.length]);
	Plot.add("circles", xpts, ypts);
	Plot.add("error bars", errbar);

	Plot.setColor("red");
	Plot.add("line", xpts, Vfit);
	
	Plot.setColor("black");
	Plot.add("line", xpts, Dfit);
	
	Plot.setColor("blue");
	Plot.add("line", xpts, Sfit);
	Plot.setLimits(0, NaN, 0, maxY*1.25);
	
	Plot.setFontSize(16, "bold");
	Plot.setColor("black");
	Plot.addText( "ID: "+d2s(j,0)+" (" +d2s(trackpts[j],0)+ " pts) :: (" +d2s(x_vals[st],0)+ "," +d2s(y_vals[st],0)+","+ d2s(slice_num[st],0)+") (x,y,z)", 0.25, 0.05);
	
	Plot.setFontSize(14, "bold");
	
//	if (VRsq>DlatRsq) {
	if (S_exp[ct] >= alpha){
		Plot.setFontSize(18, "bold");
// create a CSV list of motorized plots (corresponding sliceNo)
		if (motor) motorStr+=",";
		motor=true;
		motorStr+=d2s(sliceNo,0);
	}
	Plot.setColor("red");
	Plot.addText("Vel: "+d2s(Vel[ct],2)+umPerSec, 0.1, 0.1);
	Plot.setFontSize(12, "normal italic");
	Plot.addText("    (R="+d2s(VRsq,2)+")", 0.1, 0.13);
	Plot.setColor("black");
	Plot.setFontSize(14, "bold");
	Plot.addText("Vmax="+d2s(meanVel[j],2), 0.3, 0.16);
	
//	if (DlatRsq>VRsq) {
	if (S_exp[ct] < alpha){
		Plot.setFontSize(18, "bold");
		brownWalk=true;						// We have at least one Brownian Walk!
	}
	Plot.setColor("black");
	Plot.addText("Dlat: "+d2s(Dlat[ct],2)+um2PerSec, 0.4, 0.1);
	Plot.setFontSize(12, "normal italic");
	Plot.addText("    (R="+d2s(DlatRsq,2)+")", 0.4, 0.13);
	Plot.setColor("black");
	Plot.setFontSize(14, "bold");
	Plot.addText("Vmax="+d2s(meanVel[j],2), 0.3, 0.16);
	
	Plot.setColor("blue");
	Plot.addText("S_exp: "+d2s(S_exp[ct],2), 0.75, 0.1);
	Plot.setFontSize(12, "normal italic");	
	Plot.addText("    (R="+d2s(S_expRsq,2)+")", 0.75, 0.13);
	Plot.setColor("black");
	Plot.setFontSize(14, "bold");
	Plot.addText("Vmax="+d2s(meanVel[j],2), 0.3, 0.16);
	
	Plot.appendToStack;
	sliceNo+=1;
	
//=========== Save MSD_dT ".csv" files ===============	
	if (saveMSD){
		
		if (S_exp[ct] < alpha){
	  		for (i=0; i<minOf(range,25); i+=1){
				str= d2s(i,0) +","+ d2s(xpts[i],3) +","+ d2s(ypts[i],3)+","+ d2s(errbar[i],3);
				meanBrownianMSD[i] = meanBrownianMSD[i]+ypts[i];
				File.append(str, saveBrownian_MSD_file);
				brownCt[i]+=1;
			}

		File.append(" ", saveBrownian_MSD_file);
		}
		
		if (S_exp[ct] >= alpha){
	  		for (i=0; i<minOf(range,25); i+=1){
				str= d2s(i,0) +","+ d2s(xpts[i],3) +","+ d2s(ypts[i],3)+","+ d2s(errbar[i],3);
				meanMotorMSD[i] = meanMotorMSD[i]+ypts[i];
				File.append(str, saveMotor_MSD_file);
				motorCt[i]+=1;
			}

		File.append(" ", saveMotor_MSD_file);
		}	
	}
	
 }
 // end of "valid data" loop
 
showProgress(j, tracknd);
ct+=1;

}
 // finish the outer loop

// Now write out data to the Results log and sort the MSD vs dT plots
Plot.show();
run("Properties...", "channels=1 slices=1 frames="+d2s(nSlices(),0)+" pixel_width=1 pixel_height=1 voxel_depth=1 origin=0,0");

//split the stack into Browian Walks and Motorized trajectories
if (brownWalk) {
	rename("Bownian_Walk");
	if (motor){
		run("Make Subset...", "slices="+motorStr+" delete");
	}
}
rename("Motorized");

// blank-out unfilled columns the null rows in Results Log
for (i=nVel; i<nResults; i+=1) {
	setResult("ID_Vel", i, "");
	setResult("Vel", i, "");
	setResult("R_Vel", i, "");
	setResult("n_Vel", i, "");
}
for (i=nDlat; i<nResults; i+=1) {
	setResult("ID_Dlat", i, "");
	setResult("Dlat", i, "");
	setResult("R_Dlat", i, "");
	setResult("n_Dlat", i, "");
}
for (i=nS; i<nResults; i+=1) {
	setResult("ID_Exp", i, "");
	setResult("Exp", i, "");
	setResult("R_Exp", i, "");
	setResult("Immed_Vel", i, "");
	setResult("n_Exp", i, "");
}

updateResults();


// create the summary histograms and plots 
Plot.create("Velocity Histograms", "Velocity (um.s-1)", "Nobs");
Plot.setFontSize(16, "normal");
Plot.setLineWidth(2);
Plot.setFrameSize(256,256);
Plot.setColor("blue");

Plot.addHistogram(meanVel, 0, 0);
Plot.addText("Immed. Vel.", 0.5, 0.1);
Plot.setColor("red");
Plot.addHistogram(Vel, 0, 0);
Plot.addText("Vel.", 0.1, 0.1);
Plot.show();
run("Properties...", "channels=1 slices=1 frames=1 pixel_width=1 pixel_height=1 voxel_depth=1 origin=0,0");

Plot.create("Dlat Histogram", "Dlat (um2.s-1)", "Nobs");
Plot.setFrameSize(256,256);
Plot.setColor("blue");
Plot.addHistogram(Dlat, 0, 0);
Plot.show();
run("Properties...", "channels=1 slices=1 frames=1 pixel_width=1 pixel_height=1 voxel_depth=1 origin=0,0");

Plot.create("Exponent fit values", "Exponent", "Nobs");
Plot.setFrameSize(256,256);
Plot.setColor("blue");
Plot.addHistogram(S_exp,0.1, 0);
Plot.show();
run("Properties...", "channels=1 slices=1 frames=1 pixel_width=1 pixel_height=1 voxel_depth=1 origin=0,0");


// stick the overall meanMSD data at the end (track #999)
	if (saveMSD){
        for (i=0; i<25;i+=1){
			if (brownCt[0]>0){
				str= d2s(999,0) +","+ d2s(i*zCal,3) +","+ d2s(meanBrownianMSD[i]/brownCt[i],3);
				File.append(str, saveBrownian_MSD_file);
			}
			if (motorCt[0]>0){
				str= d2s(999,0) +","+ d2s(i*zCal,3) +","+ d2s(meanMotorMSD[i]/motorCt[i],3);
				File.append(str, saveMotor_MSD_file);
			}
		}
	}

}

//======================
function plotOverlays(){
//======================
// This is quite a large chunk of code which I put into a function 
// just to make the main program a bit more readable!
// NOTE: Arrays are "PUBLIC" so do not need to be passed to the function.

setLineWidth(lineThick); setFont("SansSerif", trackFontSize);
wM=10; wMh=wM/2; hM=10; hMh=hM/2;

st=0; nd=0;
oldtracknum=track_num[0];
showStatus("Plotting overlays");

for (j=0; j<track_num.length; j+=1){
	setColor(colors[ oldtracknum % colors.length ]);
	if (track_num[j] == oldtracknum) nd+=1;
 	else {
// paint all tracks on first frame
		i=st; 
		oldx=x_vals[i]; oldy=y_vals[i];
		while (i<nd){
			x = x_vals[i]; y = y_vals[i];
			Overlay.drawLine(oldx, oldy, x, y);
			Overlay.setPosition(1);
			oldx=x; oldy=y;
			i+=1;
		}
		st=j; nd=j;
	}
	
// paint track on each frame jumping between slices ("z") is slow but makes code simple!
// ==================================================================
	oldtracknum=track_num[j];
	i=st; 
	oldx=x_vals[i]; oldy=y_vals[i];
	while (i<nd){
		x = x_vals[i]; y = y_vals[i]; z = slice_num[j];
		Overlay.drawLine(oldx, oldy, x, y);
		Overlay.setPosition(z);
		oldx=x; oldy=y;
		i+=1;
	}

// Optional bullseye and track number label
//================================
	x = x_vals[j]; y = y_vals[j]; z = slice_num[j];
	Overlay.drawLine(oldx,oldy,x,y);
	Overlay.setPosition(z);
	if (bullseye){
		Overlay.drawLine(x-wMh,y,x-wMh+2,y);
		Overlay.setPosition(z);
		Overlay.drawLine(x+wMh,y,x+wMh-2,y);
		Overlay.setPosition(z);
		Overlay.drawLine(x,y-hMh,x,y-hMh+2);
		Overlay.setPosition(z);
		Overlay.drawLine(x,y+hMh,x,y+hMh-2);
		Overlay.setPosition(z);
		Overlay.drawEllipse(x-wMh, y-hMh, wM, hM);
		Overlay.setPosition(z);
	}
	
	if (trackID){
		Overlay.drawString(toString(track_num[j]), x-wMh, y-(hMh+2));
		Overlay.setPosition(z);
	}
	updateDisplay();
	showProgress(j,slice_num.length);
}

Overlay.add();
Overlay.show();
setSlice(1);
}

//===============================
function autoProminence(nSpots,option) {
//===============================
run("Clear Results");
print ("\\Clear");
getMinAndMax(min, max);
incr=floor(max/100);
if (incr<1) incr=1;
ct=1;
n=nSpots+1;

while ((ct<max) && (n>nSpots) && (incr>=1)) {
	do {
		run("Find Maxima...", "prominence="+d2s(ct,2)+option+" output=Count");
		n= getResult("Count", nResults()-1);
		print ("inc=",incr, "ct=",ct,"nSpots=",n);
		ct+=incr;
	} while ((ct<max) && (n>nSpots));
	ct-=incr;
	incr=incr/2;
}
run("Clear Results");
print ("\\Clear");
return (ct);
}

//=======================
function LoadMetaData() {
//=======================
//This function tries to load meta data txt file

// generate the expected metadata filename
path=getInfo("image.directory");
fName_tif=getInfo("image.filename");

errStr="Not an OME.TIF image file";
if (endsWith(fName_tif,"ome.tif")!=1) return errStr;

errStr="Ome.Tif Metadata already loaded!!";
if (!isNaN(parseFloat(Property.get("ElapsedTime-ms")))) return errStr;


fName_metadata = substring(fName_tif,0,indexOf(fName_tif,"\.ome"))+"_metadata.txt";

// See if the metadatafle exists
if (!File.exists(path+fName_metadata)) {
	errStr ="= = = = = = = = = = = = = = = = = = =\n";
	errStr+=" ========= WARNING ========= \n";
	errStr+="= = = = = = = = = = = = = = = = = = =\n";
	errStr+="Current image File: "+fName_tif+"\n";
	errStr+="No metadata file:     "+fName_metadata+"\n";
	errStr+="YOU NEED TO ENTER frame timing manually using =>Image/Property/ menu";
	return errStr;
}

// Yippee found it!.. 
// load it and parse Frame Info into all Stack Frame label fields

str = File.openAsString(path+fName_metadata);

lines=split(str, "\n");

// There is loads of "junk" metadata which we can just ignore!
// So basically just pick out the first 17 data items associated with each frame
// and store the data in the "Slice Label / Info" - then it can be accessed
// in the main program using e.g. PixSizeUm = (Property.get("PixelSizeUm");

ct=0; st=0; nd=0; frNo=0;
do {
	if ((indexOf(lines[ct], "FrameKey") > 0) && (st > 0)) nd = ct;
	if ((indexOf(lines[ct], "FrameKey") > 0) && (nd == 0)) st = ct;
	if (nd > st){
		frNo+=1;
		DoSliceLabel(st, nd, frNo,fName_tif);
		nd=0; st=0; ct-=1;
	}
	ct+=1;
} while (ct<lines.length);

// do last frame
nd = lines.length;
frNo+=1;
DoSliceLabel(st, nd, frNo,fName_tif);
		
// make best guess at average frame rate;
setSlice(1);
startT = parseFloat(Property.get("ElapsedTime-ms"))/1000;
setSlice(nSlices());
endT = parseFloat(Property.get("ElapsedTime-ms"))/1000;
interval = (endT-startT)/nSlices();

Stack.setTUnit("sec");
Stack.setFrameInterval(interval); 
Stack.setFrameRate(1/interval);

// Label each frame with exact time and scale bar
run("Remove Overlay");
setFont("SanSerif", 30, "bold");
setColor("yellow");
setSlice(1);

// get the start time
startT = parseFloat(Property.get("ElapsedTime-ms"))/1000;
for (i=0;i<nSlices();i+=1){
// ImageJ indexing starts at 1 for frames/slices
	setSlice(i+1);
	Overlay.drawString(d2s( (parseFloat(Property.get("ElapsedTime-ms"))/1000) - startT, 2 )+ " sec", 10,30);
	//t= parseInt(Property.get("Frame_t"));
	//c= parseInt(Property.get("Frame_c"));
	//z= parseInt(Property.get("Frame_z"));
	Overlay.setPosition(i+1);
	updateDisplay();
}
Overlay.add;
Overlay.show;
run("Scale Bar...", "width=50 height=0 thickness=10 font=50 color=Yellow background=None location=[Lower Right] horizontal bold overlay label");

setSlice(1);

return "Successfully Loaded metadata & set properties";
}

// ====== Set the Slice label =======
function DoSliceLabel(st,nd,frNo,fName_tif){
//		frameInfo = String.join(Array.slice(lines,st+2,nd-1),"\n");
		frameInfo = String.join(Array.slice(lines,st+1,st+19),"\n");
		frameInfo = replace(frameInfo,"\"","");
		frameInfo = replace(frameInfo,"\,","");

// Disconbobulate the "FrameKey-0-0-0: }" metadata
// Note: ImageJ indexing starts at 1 for times/colours/z's
// Micromanager uses zero indexing and tcz rather than czt.. 
		frameKey = lines[st];
		frameKey = replace(frameKey, "\"","");
		tcz = split(frameKey, "\-\:");
		frameKey = "Frame_t: "+parseInt(tcz[1])+1+"\n";
		frameKey+= "Frame_c: "+parseInt(tcz[2])+1+"\n";
		frameKey+= "Frame_z: "+parseInt(tcz[3])+1+"\n";
		frameInfo = frameKey + frameInfo + "\n";
//Add the file name as the first label
		frameInfo = fName_tif+ "\n"+ frameInfo + "\n";
  		Property.setSliceLabel(frameInfo,frNo);
}


