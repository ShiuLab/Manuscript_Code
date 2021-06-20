//java -cp ij.jar ij.ImageJ -macro macroname.ijm <image_path>
//macro "Count Seeds" 
//setBatchMode(true);
file_name = getArgument();
if (endsWith(file_name,"xls") || endsWith(file_name,"txt")){
exit();
}
//Current template
ovals = newArray(914,576,2558,2653,4048,569,2584,2679,7264,616,2568,2656,788,3994,2649,2667,4072,4016,2632,2624,7304,4010,2635,2614,872,7467,2616,2605,4104,7392,2624,2648,7296,7408,2592,2648,814,10823,2634,2624,4104,10791,2609,2632,7304,10847,2600,2607);


function count_seeds(file_name,ovals) {
        open( file_name);
		run("Set Scale...","distance=2880.0 known=60 unit=unit"); //Analyze > set scale
		setThreshold(50, 140);//Image > adjust > threshold
        path = substring(file_name,0,lastIndexOf(file_name,"\\"));
        new_path = path+'\\counts\\';
        new_base_name = substring(file_name,(lastIndexOf(file_name,'\\')+1),indexOf(file_name,'.bmp'));
        File.makeDirectory(new_path);
        savename = new_path+new_base_name+"_results";
		
        for(i=0;i<ovals.length;i+=4){
            makeOval(ovals[i+0],ovals[i+1],ovals[i+2],ovals[i+3]);
			results_name = savename+toString(i/4,0);	
			
			run("Analyze Particles...","size=0.06-infinity circularity=0.25-1.00 display clear summarize add");
			run("Flatten");
			saveAs("JPEG",savename+i/4+".jpg");
			run("Close");
			
			
			
        };
		
        saveAs("Results",savename+".xls");
        
		
        run("Close");
		
        close("*");
        exit();
}
count_seeds(file_name,ovals);

