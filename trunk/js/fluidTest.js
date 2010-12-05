var cwid;
var chei;
var imageData;
var solver;
var fluidBuffer=[];
var fw = 30;
var fh;
var lmx;
var lmy;
var isw;
var ish;
var F2C;
var C2F;
var aspectRatio2;
var frameCount = 0;
var maxInt = parseInt('FFFFFF', 32);
var canvas;
var copyContext;		/* context for "fake" canvas where we do all our pixel manipulation */
var copyCanvas;			/* fake canvas, for the main canvas to draw */
var blob;
var dragging = false;
var canpos;
var mode=1;

$(document).ready(function(e){
	canvas = document.getElementById("canvas1");	
	ctx = canvas.getContext("2d");
	cwid = parseInt(canvas.getAttribute("width"));
	chei = parseInt(canvas.getAttribute("height"));
	fh = parseInt((chei/cwid) * fw);

	copyCanvas = document.createElement("canvas");
	copyCanvas.width = fw;
	copyCanvas.height = fh;
	copyContext = copyCanvas.getContext("2d");

	F2C = cwid/(fw-1);
	C2F = (fw-5)/cwid;
	isw = 1 / cwid;
	ish = 1 / chei;

	aspectRatio2 = (cwid * ish) * (cwid * ish);	

	solver = new fluidsolver(fw, fh);	
	solver.isRgb = false;
	solver._fadeSpeed = .007;
	solver._dt = .5;
	solver._visc = .015;

	blob = new particle(cwid/2.0, cwid/2.0, 0.002, canvas, solver);
	
	imageData = copyContext.createImageData(fw, fh);
	
	canpos = $('#canvas1').position();
	lmx = canpos.left;
	lmy = canpos.top;

	$(document).mousemove(function(e) {
			handleForce(e.pageX - canpos.left, e.pageY - canpos.top);
	});
	$("#vSlider").slider({
		min: 1,
		max: 10000,
		value: solver.viscosity * 10000,
		change: function(event, ui) {
			solver.viscosity = $( "#vSlider" ).slider( "option", "value" )/10000;
		}
	});
	
	$("#fadeSlider").slider({
		min: 1000,
		max: 9990,
		value: solver.fadeSpeed * 10000,
		change: function(event, ui) {
			solver.fadeSpeed = $( "#fadeSlider" ).slider( "option", "value" )/10000;
		}
	});
	
	$("#colorSlider").slider({
		min: 0.0,
		max: 10.0,
		value: solver.colorDiffusion * 10000,
		change: function(event, ui) {
			solver.colorDiffusion = $( "#colorSlider" ).slider( "option", "value" )/10000;
		}
	});
	$("#deltaSlider").slider({
		min: 2000,
		max: 20000,
		value: solver.delta * 10000,
		change: function(event, ui) {
			solver.delta = $( "#deltaSlider" ).slider( "option", "value" )/10000;
		}
	});
	
	$( "#format" ).buttonset();
	$("#wrapx").click(function() {
		solver.wrap_x = !solver.wrap_x;
	});
	$("#wrapy").click(function() {
		solver.wrap_y = !solver.wrap_y;
	});
	
	$( "#modeButton" ).button();
	
	$('#modeButton').mousedown(function(e) {
		mode ++;
		if (mode>1) mode=0;
		ctx.fillStyle="#000000"
		ctx.fillRect(0, 0, cwid, chei);
		ctx.fill();
	});
	
	var timer = setInterval ( blit, 40 );
});

function mouseUpHandler(e) { dragging = false; }


function handleForce(x, y) {
	addForce(x * isw, y * ish, (x - lmx) * isw, (y - lmy) * ish);
	lmx = x;
	lmy = y;
}

function blit() {
	solver.update();
	if (mode==0) drawCircles();
	else if (mode==1) {
		drawFluidBitmap();
		frameCount = ++frameCount % maxInt;
	}
}

function drawFluidBitmap() {
	var fi;
	var li;
	var index = 0;
	
	var dat = imageData.data;
	
	for(var j = 1; j < fh-1; ++j) {
		for(var i = 1; i < fw-1; ++i) {
			fi = parseInt(i + fw * j);
			li = fi*4;
			dat[li]= (solver.r[fi] * maxInt/16);
			dat[li+1]= (solver.g[fi] * maxInt/16);
			dat[li+2]= (solver.b[fi] * maxInt/16);
			dat[li+3]=maxInt;
		}
	}
	
	
	imageData.data = dat;
	copyContext.putImageData(imageData, 0, 0);
	ctx.drawImage(copyCanvas, 0, 0, cwid, chei);
}

function drawCircles() {
	ctx.fillStyle="#000000"
	ctx.fillRect(0, 0, cwid, chei);
	ctx.beginPath();
	var fi;
	var li;
	var rad = 0;
	for(var j = 0; j < fh; ++j) {
		for(var i = 0; i < fw; ++i) {
			fi = parseInt(i + fw * j);
			rad = (solver.r[fi]+solver.g[fi]+solver.b[fi])*F2C;
			ctx.moveTo(i * F2C + rad, j * F2C);
			ctx.arc(i * F2C, j * F2C, rad, 0, Math.PI * 2, false);
		}
	}
	ctx.closePath();
	ctx.strokeStyle = "#73C7FF";
	ctx.lineWidth = 3;
	ctx.fillStyle = "#43647A";
	ctx.stroke();
	ctx.fill();
	//drawBlob();
}
/*
function drawBlob() {
	var rad = 20;
	
	ctx.beginPath();
	var tx = blob.x * C2F;
	var ty = blob.y * C2F;
	var fi = parseInt(tx + fw * ty);
	blob.addForce(solver.u[fi], solver.v[fi]);
	blob.update();
	ctx.moveTo(blob.x + rad, blob.y);
	ctx.arc(blob.x, blob.y, rad, 0, Math.PI * 2, false);
	ctx.closePath();
	ctx.strokeStyle = "#FF0000";
	ctx.lineWidth = 3;
	ctx.fillStyle = "#FF6666";
	ctx.stroke();
	ctx.fill();	
}
*/
function addForce(x, y, dx, dy) {
	var speed = dx * dx  + dy * dy * aspectRatio2;

	if(speed > 0) {
		if (x < 0) x = 0;
		else if (x > 1) x = 1;
		if (y < 0) y = 0;
		else if (y > 1) y = 1;

		var index = solver.getIndexForNormalizedPosition(x, y);
		
		
		var colorMult  = 30;
		var velocityMult  = 20.0;

		var hue = ((x + y) * 180 + frameCount) % 360;
		var rgb = HSB2GRB(hue, 1, 1);

		solver.rOld[index]  += rgb.r * colorMult;
		solver.gOld[index]  += rgb.g * colorMult;
		solver.bOld[index]  += rgb.b * colorMult;
		
		solver.uOld[index] += dx * velocityMult;
		solver.vOld[index] += dy * velocityMult;
	}
}

function HSB2GRB(h, s, b) {	
	h = parseInt(h) % 360;
	var i = parseInt(parseInt(h / 60.0) % 6);
	var f = h / 60.0 - parseInt(h / 60.0);
	var p = b * (1 - s);
	var q = b * (1 - s * f);
	var t = b * (1 - (1 - f) * s);
	switch (i) {   
		case 0: return { r:b, g:t, b:p };
		case 1: return { r:q, g:b, b:p };
		case 2: return { r:p, g:b, b:t };
		case 3: return { r:p, g:q, b:b };
		case 4: return { r:t, g:p, b:b };
		case 5: return { r:b, g:p, b:q };
	}
	return { r:0, g:0, b:0 };
}