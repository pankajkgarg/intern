
var settings = {
	"width": 1200,
	"height": 700,
	"padding" : [1, 1, 1, 1], // top right bottom left
	"fontResize": false,
	"circleResize": false,
	"displayZoom": 1.8,
};

var xTemp,yTemp, zTemp;


//preprocessing data
for(var i = 0; i < data.length; i++) {
	data[i].x = Number(data[i].x);
	data[i].y = Number(data[i].y);
	
}

//Finding max and min value of x and y
var minX = 10000, minY = 10000, maxX = -100000, maxY = -100000;
var minGenes = 100000, maxGenes = 0;
var minScore = 10000, maxScore = -10000;
for(var i = 0; i < data.length; i++) {
	if (data[i].x > maxX) { maxX = data[i].x };
	if (data[i].x < minX) { minX = data[i].x };
	if (data[i].y > maxY) maxY = data[i].y;
	if (data[i].y < minY) minY = data[i].y;
	
	if (data[i].score > maxScore) maxScore = data[i].score;
	if (data[i].score < minScore) minScore = data[i].score;
	if (data[i].numGenes > maxGenes) maxGenes = data[i].numGenes;
	if (data[i].numGenes < minGenes) minGenes = data[i].numGenes;
}


var x = d3.scale.linear()
	.domain([minX * 1.7, maxX * 1.5])
	.range([0, settings.width]);
	
// y scale - inverted domain
var y = d3.scale.linear()
	.domain([maxY * 1.7,  minY * 1.5])
	.range([0, settings.height]);


for(var i = 0; i < data.length; i++) {
	data[i]["r"] = x(data[i].score);
}


tx = function(d) { return "translate(" + x(d) + ",0)"; },
ty = function(d) { return "translate(0," + y(d) + ")"; },
stroke = function(d) { return d ? "#ccc" : "#666"; };

// var chart = d3.select("body")
	// .append("svg:svg")
		// .attr("class", "chart")
		// .attr("width", settings.width)
		// .attr("height", settings.height);
		


		

//color range
var color = d3.scale.linear()
	.domain([minScore, maxScore])
	.range(["white", "darkblue"]);
	
var sizeScale = d3.scale.log()
	.domain([minGenes, maxGenes])
	.range([3,12]);
	
//color = d3.scale.category10();
var symbol = d3.scale.ordinal().range(d3.svg.symbolTypes);

var baseSVG = d3.select("#svgRoot")
	.append("svg:svg")
		.attr("width", settings.width + settings.padding[1] + settings.padding[3])
		.attr("height", settings.height + settings.padding[0] + settings.padding[2])
		.call(d3.behavior.zoom()
		.on("zoom", redraw))


		
frameworkSVG = baseSVG.append("svg:g")
		.attr("render-order", 10)
		.attr("transform", "translate(" + settings.padding[3] + "," + settings.padding[0] + ")")
		
frameworkSVG.append("svg:rect")
		.attr("width", settings.width)
		.attr("height", settings.height)
		.attr("stroke", stroke)
		.attr("fill", "none");

svg = baseSVG.append("svg:g")
		.attr("render-order", 0)
		.attr("transform", "translate(" + settings.padding[3] + "," + settings.padding[0] + ")")
	
radiusFunc = function(d) {
	return sizeScale(d.numGenes);//d.r * 0.005;
}	
		
svg.selectAll("circle")
	.data(data)
.enter().append("svg:circle")
	.attr("class", "dot")
	.attr("stroke", function(d, i) {return color(d.score) })
	.attr("fill", function(d, i) { return color(d.score)})
	.attr("fill-opacity", 0.2)
	.attr("cx", function(d, i) { return x(d.x); })
	.attr("cy", function(d, i) { return y(d.y); })
	.attr("r", radiusFunc);
	;
	
svg.selectAll("text")
	.data(data)
	.enter().append("svg:text")
		.attr("class", "label")
		.attr("x", function(d) { return x(d.x); }) 	//  + (x(d.z) * 0.001)
		.attr("y", function(d) { return y(d.y); })
		.attr("text-anchor", "middle")
		//.text(function(d) { return String(d.label);});



maxLevel = 0
var allLabels = [];
for(var level in zoomLabels){
	var labels = zoomLabels[level]
	for (var i = 0; i < labels.length; i++){
		labels[i]["zoom"] = level;
		labels[i]["r"] = 20000/level;
		allLabels.push(labels[i])
	}
	if (Number(level) > Number(maxLevel)) {	
		maxLevel = Number(level)
	}
}

var maxLabelScore = 0, minLabelScore = 100000;
for (var i = 0; i < allLabels.length; i++){
	if (allLabels[i].score > maxLabelScore)  maxLabelScore = allLabels[i].score
	if (allLabels[i].score < minLabelScore)  minLabelScore = allLabels[i].score
}
var labelSize = d3.scale.linear()
	.domain([minLabelScore, maxLabelScore])
	.range([14, 45])

svg.selectAll("text.zoomLabel")
	.data(allLabels)
	.enter().append("svg:text")
	.attr("class", function(d) { return "zoomLabel invisible level_" + String(d.zoom) })
	.attr("x", function(d) { return x(d.x); })
	.attr("y", function(d) { return y(d.y); })
	.attr("font-size", function(d) { return labelSize(d.score);  })
	.attr("text-anchor", "middle")
	.text(function(d){ return String(d.label); })



lastValue = false
function displayZoomLabels(level) {
	$( "#amount" ).val( level);
	if (level != lastValue){
		if (lastValue != false){
			d3.selectAll("text.visible").attr("class", "zoomLabel invisible level_" + String(lastValue) )
		}
		d3.selectAll("text.level_" + String(level)).attr("class", "zoomLabel visible")
		lastValue = level
	}
}

$(function() {
	$( "#slider" ).slider({
		value:1,
		min: 1,
		max: maxLevel,
		step: 1,
		slide: function( event, ui ) {
			displayZoomLabels(ui.value);
			
			svg.selectAll("text.zoomLabel")
				.attr("x", function(d) {  return x(d.x);})		
				.attr("y", function(d) {return y(d.y); })		
		}
	});
	$( "#amount" ).val( $( "#slider" ).slider( "value" ) );
	
});
displayZoomLabels(1)

var zoomScale = d3.scale.linear()
	.domain([1.0, 3.0])
	.range([1, Number(maxLevel)])
	
	
redraw();

var lastScale = false;
function redraw() {
	if (d3.event) d3.event.transform(x, y);

	var fx = x.tickFormat(10),
	  fy = y.tickFormat(10);

	// Regenerate x-ticks…
	var gx = frameworkSVG.selectAll("g.x")
	  .data(x.ticks(10), String)
	  .attr("transform", tx)
	  .attr("class", "x");

	gx.select("text")
	  .text(fx);

	var gxe = gx.enter().insert("svg:g", "rect")
	  .attr("class", "x")
	  .attr("transform", tx);

	gxe.append("svg:line")
	  .attr("stroke", stroke)
	  .attr("y1", 0)
	  .attr("y2", settings.height);

	gxe.append("svg:text")
	  .attr("y", settings.height)
	  .attr("dy", "1em")
	  .attr("text-anchor", "middle")
	  .text(fx);

	gx.exit().remove();

	// Regenerate y-ticks…
	var gy = frameworkSVG.selectAll("g.y")
	  .data(y.ticks(10), String)
	  .attr("transform", ty)
	  .attr("class", "y");

	gy.select("text")
	  .text(fy);

	var gye = gy.enter().insert("svg:g", "rect")
	  .attr("class", "y")
	  .attr("transform", ty);

	gye.append("svg:line")
	  .attr("stroke", stroke)
	  .attr("x1", 0)
	  .attr("x2", settings.width);

	gye.append("svg:text")
	  .attr("x", -3)
	  .attr("dy", ".35em")
	  .attr("text-anchor", "end")
	  .text(fy);

	gy.exit().remove();
	
	var s = x(1) - x(0)
	var scale = (d3.event ? d3.event.scale : 1);
	// This equation give zoom levels as 1, 2, 3 instead of 2, 4, 8, 16
	var zoomLevel = (d3.event ? (Math.log(d3.event.scale)/ Math.LN2) : 1) 
	
	//console.log(zoomLevel);
	//console.log("scale ", scale)
	
	//Display zoom labels corresponding to current scale
	if (scale == 0) {
		scale = 1;
	}
	
	scale = Math.ceil(zoomScale(zoomLevel))
	displayZoomLabels(scale)
	//console.log("Transformed scale ", scale)
	
	
//	if (scale != lastScale){
//		if (lastScale != false){
//			d3.selectAll("text.visible").attr("class", "zoomLabel invisible level_" + String(lastScale) )
//		}
//		d3.selectAll("text.level_" + String(scale)).attr("class", "zoomLabel visible")
//		lastScale = scale
//	}
	svg.selectAll("text.zoomLabel")
		.attr("x", function(d) {  return x(d.x); })		
		.attr("y", function(d) {return y(d.y); })		
	
	
	var allowedLabel = 80/zoomLevel
	
	svg.selectAll("circle")
		//.attr("cx", function(d, i) { tempX = x(d.x); return ( (tempX > 1) ? tempX :  "-200"); })
		//.attr("cy", function(d, i) { tempY = y(d.y); return ( (tempY < settings.height) ? tempY : "-200"); })
		.attr("cx", function(d) { return x(d.x) ;  })
		.attr("cy", function(d) { return y(d.y); })
		
		;
	
	

	
	
	
	svg.selectAll("text.label")
		.attr("x", function(d) {  return x(d.x) + radiusFunc(d) ; })		//  + (x(d.z) * 0.001) + 2
		.attr("y", function(d) {return y(d.y) - radiusFunc(d) - 5 ; })		//  + (0.5 * x(d.z) * 0.001)
		.text(function(d) { return ((zoomLevel >= settings.displayZoom)? d.label : "") ;});  
	
	if (settings.fontResize) { svg.selectAll("text.label").attr("font-size", function(d) { return s * 0.003; }); }
	if (settings.circleResize) { svg.selectAll("circle").attr("r", function(d) { tempX =  s * d.score * 0.002;   return tempX;}) ; }

		
}

