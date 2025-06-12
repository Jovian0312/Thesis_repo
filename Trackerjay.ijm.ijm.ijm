//Clear previous results
run("Clear Results");

//Get total number of frames
getDimensions(width, height, channels, slices, frames);

//Loop through each frame
for(i = 1; i <= frames; i++) {
	setSlice(i);//Move to the i-th frame
	run("Plot Profile");
	
	//Get number of points in profile
	n = nResults;
	
	//Initialize max values
	maxVal = -1;
	maxPos = -1;
	
	//Search for max intensity and its position
	for (j = 0; j < n; j++) {
		val = getResults("Y", j);
		if( val > maxVal ) {
			maxVal = val;
			maxPos = getResult("X", j);
		
		}
	}
	
	
	// Store current results
	setResult("Frame", i-1, i);
	setResult("Max Position", i-1, maxPos);
	setResult("Max Intensity", i-1, maxVal);
	updateResults();
	

}
saveAs("Results","~/Desktop/bead_tracking_result.csv")