
function fluidsolver(width, height) {

	/***********************************************************************
 
	* This is a class for solving real-time fluid dynamics simulations based on Navier-Stokes equations
	* and ported from Memo Akten's AS3 vertsion at http://code.google.com/p/in-spirit/source/browse/#svn/trunk/projects/FluidSolver
	
	***********************************************************************/

	this.FLUID_DEFAULT_NX = 50;
	this.FLUID_DEFAULT_NY = 50;
	this.FLUID_DEFAULT_DT = 1.0;
	this.FLUID_DEFAULT_VISC = 0.0001;
	this.FLUID_DEFAULT_COLOR_DIFFUSION = 0.0;
	this.FLUID_DEFAULT_FADESPEED = 0.3;
	this.FLUID_DEFAULT_SOLVER_ITERATIONS = 10;
	this.FLUID_DEFAULT_VORTICITY_CONFINEMENT = false;

	this.delta = this.FLUID_DEFAULT_DT;
	this.fadeSpeed = this.FLUID_DEFAULT_FADESPEED;
	var _solverIterations = this.FLUID_DEFAULT_SOLVER_ITERATIONS;
	this.colorDiffusion = this.FLUID_DEFAULT_COLOR_DIFFUSION;
	var _doVorticityConfinement = this.FLUID_DEFAULT_VORTICITY_CONFINEMENT;


	this.r = [];
	this.g = [];
	this.b = [];
	this.u = [];
	this.v = [];
	this.rOld = [];
	this.gOld = [];
	this.bOld = [];
	this.uOld = [];
	this.vOld = [];

	var curl_abs = [];
	var curl_orig = [];

	var _NX;
	var _NY;
	var _NX2;
	var _NY2;
	var _invNumCells;

	this.isRGB = true;				// for monochrome, only update r

	this.wrap_x = false;
	this.wrap_y = false;
	this.viscosity = 0.0;	
	var _tmp = [];
	var _avgDensity = 0;			// this will hold the average color of the last frame (how full it is)
	var _uniformity = 0;			// this will hold the uniformity of the last frame (how uniform the color is);
	var _avgSpeed = 0;
	
	_NX = width-2;
	_NY = height;
	_NX2 = _NX + 2;
	_NY2 = _NY + 2;
	this.numCells = _NX2 * _NY2;
	_invNumCells = 1.0 / this.numCells;
	
	this.width = _NX2;
	this.height = _NY2;
	
	this.reset = function() {
		var fixed = false;
		this.r = [];
		this.g = [];
		this.b = [];
		this.u = [];
		this.v = [];
		this.rOld = [];
		this.gOld = [];
		this.bOld = [];
		this.uOld = [];
		this.vOld = [];
		curl_abs = [];
		curl_orig = [];
		
		var i = this.numCells;
		while ( --i > -1 ) {
			this.u[i] = this.uOld[i] = this.v[i] = this.vOld[i] = 0.0;
			this.r[i] = this.rOld[i] = this.g[i] = this.gOld[i] = this.b[i] = this.bOld[i] = 0;
			curl_abs[i] = curl_orig[i] = 0;
		}		
	}
	this.reset();

	this.update = function() {
		this.addSourceUV();

		if( _doVorticityConfinement )
		{
			this.calcVorticityConfinement(this.uOld, this.vOld);
			this.addSourceUV();
		}

		this.swapUV();
		this.diffuseUV(this.viscosity);

		this.project(this.u, this.v, this.uOld, this.vOld);

		this.swapUV();

		this.advect(1, this.u, this.uOld, this.uOld, this.vOld);
		this.advect(2, this.v, this.vOld, this.uOld, this.vOld);

		this.project(this.u, this.v, this.uOld, this.vOld);

		if(this.isRGB) {
			this.addSourceRGB();
			this.swapRGB();

			if( this.colorDiffusion != 0 && this.delta != 0 )
	        {
				this.diffuseRGB(this.colorDiffusion);
				this.swapRGB();
	        }

			this.advectRGB(this.u, this.v);

			this.fadeRGB();
		} else {
			this.addSource(this.r, this.rOld);
			this.swapR();

			if( this.colorDiffusion != 0 && this.delta != 0 )
	        {
				this.diffuse(0, this.r, this.rOld, this.colorDiffusion);
				this.swapRGB();
	        }

			this.advect(0, this.r, this.rOld, this.u, this.v);	
			this.fadeR();
		}		
	}

	this.calcVorticityConfinement = function(_x, _y) {
		var i;
		var j;
		var dw_dx;
		var dw_dy;
		var length;
		var index;
		var vv;

		for (j = _NY; j > 0; --j) {
			index = FLUID_IX(_NX, j);
			for (i = _NX; i > 0; --i)
			{
				dw_dy = this.u[parseInt(index + _NX2)] - this.u[parseInt(index - _NX2)];
				dw_dx = this.v[parseInt(index + 1)] - this.v[parseInt(index - 1)];

				vv = (dw_dy - dw_dx) * .5;

				curl_orig[ index ] = vv;
				curl_abs[ index ] = vv < 0 ? -vv : vv;

				--index;
			}
		}

		for (j = _NY-1; j > 1; --j)
		{
			index = FLUID_IX(_NX-1, j);
			for (i = _NX-1; i > 1; --i)
			{
				dw_dx = curl_abs[parseInt(index + 1)] - curl_abs[parseInt(index - 1)];
				dw_dy = curl_abs[parseInt(index + _NX2)] - curl_abs[parseInt(index - _NX2)];

				length = Math.sqrt(dw_dx * dw_dx + dw_dy * dw_dy) + 0.000001;

				length = 2 / length;
				dw_dx *= length;
				dw_dy *= length;

				vv = curl_orig[ index ];

				_x[ index ] = dw_dy * -vv;
				_y[ index ] = dw_dx * vv;

				--index;
			}
		}
	}

	this.fadeR = function() {
		var holdAmount = 1 - this.fadeSpeed;

		_avgDensity = 0;
		_avgSpeed = 0;

		var totalDeviations = 0;
		var currentDeviation;
		var tmp_r;

		var i = this.numCells;
		while ( --i > -1 ) {
			// clear old values
			this.uOld[i] = this.vOld[i] = 0; 
			this.rOld[i] = 0;

			// calc avg speed
			_avgSpeed += this.u[i] * this.u[i] + this.v[i] * this.v[i];

			// calc avg density
			tmp_r = Math.min(1.0, this.r[i]);
			_avgDensity += tmp_r;	// add it up

			// calc deviation (for uniformity)
			currentDeviation = tmp_r - _avgDensity;
			totalDeviations += currentDeviation * currentDeviation;

			// fade out old
			this.r[i] = tmp_r * holdAmount;
		}
		_avgDensity *= _invNumCells;

		_uniformity = 1.0 / (1 + totalDeviations * _invNumCells);		// 0: very wide distribution, 1: very uniform		
	}

	this.fadeRGB = function() {
		var holdAmount = 1 - this.fadeSpeed;

		_avgDensity = 0;
		_avgSpeed = 0;

		var totalDeviations = 0;
		var currentDeviation;
		var density;

		var tmp_r;
		var tmp_g;
		var tmp_b;

		var i = this.numCells;
		while ( --i > -1 ) {
			// clear old values
			this.uOld[i] = this.vOld[i] = 0; 
			this.rOld[i] = 0;
			this.gOld[i] = this.bOld[i] = 0;

			// calc avg speed
			_avgSpeed += this.u[i] * this.u[i] + this.v[i] * this.v[i];

			// calc avg density
			tmp_r = Math.min(1.0, this.r[i]);
			tmp_g = Math.min(1.0, this.g[i]);
			tmp_b = Math.min(1.0, this.b[i]);

			density = Math.max(tmp_r, Math.max(tmp_g, tmp_b));
			_avgDensity += density;	// add it up

			// calc deviation (for uniformity)
			currentDeviation = density - _avgDensity;
			totalDeviations += currentDeviation * currentDeviation;

			// fade out old
			this.r[i] = tmp_r * holdAmount;
			this.g[i] = tmp_g * holdAmount;
			this.b[i] = tmp_b * holdAmount;

		}
		_avgDensity *= _invNumCells;
		_avgSpeed *= _invNumCells;

		_uniformity = 1.0 / (1 + totalDeviations * _invNumCells);		// 0: very wide distribution, 1: very uniform		
	}

	this.addSourceUV = function () {
		var i = this.numCells;
		while ( --i > -1 ) {
			this.u[i] += this.delta * this.uOld[i];
			this.v[i] += this.delta * this.vOld[i];
		}
	}

	this.addSourceRGB = function() {
		var i = this.numCells;
		while ( --i > -1 ) {
			this.r[i] += this.delta * this.rOld[i];
			this.g[i] += this.delta * this.gOld[i];
			this.b[i] += this.delta * this.bOld[i];		
		}	
	}

	this.addSource = function(x, x0) {
		var i = this.numCells;
		while ( --i > -1 ) {
			x[i] += this.delta * x0[i];
		}
	}

	this.advect = function(b, _d, d0, du, dv) {
		var i;
		var j;
		var i0;
		var j0;
		var i1;
		var j1;
		var index;
		var x;
		var y;
		var s0;
		var t0;
		var s1;
		var t1;
		var dt0x = this.delta * _NX;
		var dt0y = this.delta * _NY;

		for (j = _NY; j > 0; --j) {
			for (i = _NX; i > 0; --i) {

				index = FLUID_IX(i, j);

				x = i - dt0x * du[index];
				y = j - dt0y * dv[index];

				if (x > _NX + 0.5) x = _NX + 0.5;
				if (x < 0.5) x = 0.5;

				i0 = parseInt(x);
				i1 = i0 + 1;

				if (y > _NY + 0.5) y = _NY + 0.5;
				if (y < 0.5) y = 0.5;

				j0 = parseInt(y);
				j1 = j0 + 1;

				s1 = x - i0;
				s0 = 1 - s1;
				t1 = y - j0;
				t0 = 1 - t1;

				_d[index] = s0 * (t0 * d0[FLUID_IX(i0, j0)] + t1 * d0[FLUID_IX(i0, j1)]) + s1 * (t0 * d0[FLUID_IX(i1, j0)] + t1 * d0[FLUID_IX(i1, j1)]);

			}
		}
		this.setBoundary(b, _d);		
	}

	this.advectRGB = function(du, dv) {
		var i;
		var j;
		var i0;
		var j0;
		var x;
		var y;
		var s0;
		var t0;
		var s1;
		var t1;
		var index;
		var dt0x = this.delta * _NX;
		var dt0y = this.delta * _NY;

		for (j = _NY; j > 0; --j) 
		{
			for (i = _NX; i > 0; --i)
			{
				index = FLUID_IX(i, j);
				x = i - dt0x * du[index];
				y = j - dt0y * dv[index];

				if (x > _NX + 0.5) x = _NX + 0.5;
				if (x < 0.5)     x = 0.5;

				i0 = parseInt(x);

				if (y > _NY + 0.5) y = _NY + 0.5;
				if (y < 0.5)     y = 0.5;

				j0 = parseInt(y);

				s1 = x - i0;
				s0 = 1 - s1;
				t1 = y - j0;
				t0 = 1 - t1;


				i0 = FLUID_IX(i0, j0);
	            j0 = i0 + _NX2;
	            this.r[index] = s0 * ( t0 * this.rOld[i0] + t1 * this.rOld[j0] ) + s1 * ( t0 * this.rOld[parseInt(i0+1)] + t1 * this.rOld[parseInt(j0+1)] );
	            this.g[index] = s0 * ( t0 * this.gOld[i0] + t1 * this.gOld[j0] ) + s1 * ( t0 * this.gOld[parseInt(i0+1)] + t1 * this.gOld[parseInt(j0+1)] );                  
	            this.b[index] = s0 * ( t0 * this.bOld[i0] + t1 * this.bOld[j0] ) + s1 * ( t0 * this.bOld[parseInt(i0+1)] + t1 * this.bOld[parseInt(j0+1)] );				
			}
		}
		this.setBoundaryRGB();		
	}

	this.diffuse = function(b, c, c0, _diff) {
		var a = this.delta * _diff * _NX * _NY;
		this.linearSolver(b, c, c0, a, 1.0 + 4 * a);
	}

	this.diffuseRGB = function(_diff) {
		var a = this.delta * _diff * _NX * _NY;
		this.linearSolverRGB(a, 1.0 + 4 * a);
	}

	this.diffuseUV = function(_diff) {
		var a = this.delta * _diff * _NX * _NY;
		this.linearSolverUV(a, 1.0 + 4 * a);
	}

	this.project = function(x, y, p, div) {
		var i;
		var j;
		var index;

		var h = -0.5 / _NX;

		for (j = _NY; j > 0; --j) 
	    {
			index = FLUID_IX(_NX, j);
			for (i = _NX; i > 0; --i)
			{
				div[index] = h * ( x[parseInt(index+1)] - x[parseInt(index-1)] + y[parseInt(index+_NX2)] - y[parseInt(index-_NX2)] );
				p[index] = 0;
				--index;
			}
	    }

		this.setBoundary(0, div);
		this.setBoundary(0, p);

		this.linearSolver(0, p, div, 1, 4);

		var fx = 0.5 * _NX;
		var fy = 0.5 * _NY;
		for (j = _NY; j > 0; --j) 
		{
			index = FLUID_IX(_NX, j);
			for (i = _NX; i > 0; --i)
			{
				x[index] -= fx * (p[parseInt(index+1)] - p[parseInt(index-1)]);
				y[index] -= fy * (p[parseInt(index+_NX2)] - p[parseInt(index-_NX2)]);
				--index;
			}
		}

		this.setBoundary(1, x);
		this.setBoundary(2, y);

	}

	this.linearSolver = function(b, x, x0, a, c) {
		var k;
		var i;
		var j;
		var index;

		if( a == 1 && c == 4 )
		{
			for (k = 0; k < _solverIterations; ++k) 
			{
				for (j = _NY; j > 0 ; --j) 
				{
					index = FLUID_IX(_NX, j);
					for (i = _NX; i > 0 ; --i)
					{
						x[index] = ( x[parseInt(index-1)] + x[parseInt(index+1)] + x[parseInt(index - _NX2)] + x[parseInt(index + _NX2)] + x0[index] ) * .25;
						--index;                                
					}
				}
				this.setBoundary( b, x );
			}
		}
		else
		{
			c = 1 / c;
			for (k = 0; k < _solverIterations; ++k) 
			{
				for (j = _NY; j > 0 ; --j) 
				{
					index = FLUID_IX(_NX, j);
					for (i = _NX; i > 0 ; --i)
					{
						x[index] = ( ( x[parseInt(index-1)] + x[parseInt(index+1)] + x[parseInt(index - _NX2)] + x[parseInt(index + _NX2)] ) * a + x0[index] ) * c;
						--index;
					}
				}
				this.setBoundary( b, x );
			}
		}		
	}

	this.linearSolverRGB = function(a, c) {
		var k;
		var i;
		var j;
		var index3;
		var index4;
		var index;
		c = 1 / c;
		for ( k = 0; k < _solverIterations; ++k )
		{           
		    for (j = _NY; j > 0; --j)
		    {
		            index = FLUID_IX(_NX, j );
					index3 = index - _NX2;
					index4 = index + _NX2;
					for (i = _NX; i > 0; --i)
					{       
						this.r[index] = ( ( this.r[parseInt(index-1)] + this.r[parseInt(index+1)]  +  this.r[index3] + this.r[index4] ) * a  +  this.rOld[index] ) * c;
						this.g[index] = ( ( this.g[parseInt(index-1)] + this.g[parseInt(index+1)]  +  this.g[index3] + this.g[index4] ) * a  +  this.gOld[index] ) * c;
						this.b[index] = ( ( this.b[parseInt(index-1)] + this.b[parseInt(index+1)]  +  this.b[index3] + this.b[index4] ) * a  +  this.bOld[index] ) * c;                                

						--index;
						--index3;
						--index4;
					}
			}
			this.setBoundaryRGB();
		}		
	}

	this.linearSolverUV = function(a, c) {
		var index;
		var k;
		var i;
		var j;
		c = 1 / c;
		for (k = 0; k < _solverIterations; ++k) {
			for (j = _NY; j > 0; --j) {
				index = FLUID_IX(_NX, j);
				for (i = _NX; i > 0; --i) {
					this.u[index] = ( ( this.u[parseInt(index-1)] + this.u[parseInt(index+1)] + this.u[parseInt(index - _NX2)] + this.u[parseInt(index + _NX2)] ) * a  +  this.uOld[index] ) * c;
					this.v[index] = ( ( this.v[parseInt(index-1)] + this.v[parseInt(index+1)] + this.v[parseInt(index - _NX2)] + this.v[parseInt(index + _NX2)] ) * a  +  this.vOld[index] ) * c;
					--index;
				}
			}
			this.setBoundary( 1, this.u );
	        this.setBoundary( 2, this.v );
		}

	}

	this.setBoundary = function(bound, x) {
		var dst1 = FLUID_IX(0, 1);
		var dst2 = FLUID_IX(_NX+1, 1 );
		var src1 = FLUID_IX(1, 1);
		var src2 = FLUID_IX(_NX, 1);
		var i;
		var step = FLUID_IX(0, 1) - FLUID_IX(0, 0);

		if( this.wrap_x ) {
			src1 ^= src2;
			src2 ^= src1;
			src1 ^= src2;
		}
		if( bound == 1 && !this.wrap_x ) {
			for (i = _NY; i > 0; --i )
			{
				x[dst1] = -x[src1];     dst1 += step;   src1 += step;   
				x[dst2] = -x[src2];     dst2 += step;   src2 += step;   
			}
		} else {
			for (i = _NY; i > 0; --i )
			{
				x[dst1] = x[src1];      dst1 += step;   src1 += step;   
				x[dst2] = x[src2];      dst2 += step;   src2 += step;   
			}
		}

		dst1 = FLUID_IX(1, 0);
		src1 = FLUID_IX(1, 1);
		dst2 = FLUID_IX(1, _NY+1);
		src2 = FLUID_IX(1, _NY);

		if( this.wrap_y ) {
			src1 ^= src2;
			src2 ^= src1;
			src1 ^= src2;
		}
		if( bound == 2 && !this.wrap_y ) {
			for (i = _NX; i > 0; --i )
			{
			        x[dst1++] = -x[src1++]; 
			        x[dst2++] = -x[src2++]; 
			}
		} else {
			for (i = _NX; i > 0; --i )
			{
			        x[dst1++] = x[src1++];
			        x[dst2++] = x[src2++];  
			}
		}

		x[FLUID_IX(  0,   0)] = 0.5 * (x[FLUID_IX(1, 0  )] + x[FLUID_IX(  0, 1)]);
		x[FLUID_IX(  0, _NY+1)] = 0.5 * (x[FLUID_IX(1, _NY+1)] + x[FLUID_IX(  0, _NY)]);
		x[FLUID_IX(_NX+1,   0)] = 0.5 * (x[FLUID_IX(_NX, 0  )] + x[FLUID_IX(_NX+1, 1)]);
		x[FLUID_IX(_NX+1, _NY+1)] = 0.5 * (x[FLUID_IX(_NX, _NY+1)] + x[FLUID_IX(_NX+1, _NY)]);		
	}

	this.setBoundaryRGB = function() {
		if( !this.wrap_x && !this.wrap_y ) return;
		var dst1;
		var dst2;
		var src1;
		var src2;
		var i;
		var step = FLUID_IX(0, 1) - FLUID_IX(0, 0);

		if ( this.wrap_x ) {
			dst1 = FLUID_IX(0, 1);
			src1 = FLUID_IX(1, 1);
			dst2 = FLUID_IX(_NX+1, 1 );
			src2 = FLUID_IX(_NX, 1);

			src1 ^= src2;
			src2 ^= src1;
			src1 ^= src2;

			for (i = _NY; i > 0; --i )
			{
				this.r[dst1] = this.r[src1]; this.g[dst1] = this.g[src1]; this.b[dst1] = this.b[src1]; dst1 += step;   src1 += step;   
				this.r[dst2] = this.r[src2]; this.g[dst2] = this.g[src2]; this.b[dst2] = this.b[src2]; dst2 += step;   src2 += step;   
			}
		}

		if ( this.wrap_y ) {
			dst1 = FLUID_IX(1, 0);
			src1 = FLUID_IX(1, 1);
			dst2 = FLUID_IX(1, _NY+1);
			src2 = FLUID_IX(1, _NY);

			src1 ^= src2;
			src2 ^= src1;
			src1 ^= src2;

			for (i = _NX; i > 0; --i )
			{
				this.r[dst1] = this.r[src1]; this.g[dst1] = this.g[src1]; this.b[dst1] = this.b[src1];  ++dst1; ++src1;   
				this.r[dst2] = this.r[src2]; this.g[dst2] = this.g[src2]; this.b[dst2] = this.b[src2];  ++dst2; ++src2;   
			}
		}		
	}

	this.swapUV = function() {
		_tmp = this.u; 
		this.u = this.uOld; 
		this.uOld = _tmp;

		_tmp = this.v; 
		this.v = this.vOld; 
		this.vOld = _tmp;
	}

	this.swapR = function() { 
		_tmp = this.r;
		this.r = this.rOld;
		this.rOld = _tmp;
	}

	this.swapRGB = function() { 
		_tmp = this.r;
		this.r = this.rOld;
		this.rOld = _tmp;

		_tmp = this.g;
		this.g = this.gOld;
		this.gOld = _tmp;

		_tmp = this.b;
		this.b = this.bOld;
		this.bOld = _tmp;
	}

	function FLUID_IX(i, j) { 
		return parseInt(i + _NX2 * j);
	}

	/**
	*	IMPLEMENT THESE LATER, IF NECCESARY
	**/

	/*
	function shiftLeft() {

	}

	function shiftRight() {

	}

	function shiftUp() {

	}

	function shiftDown() {

	}
	*/

	this.randomizeColor = function() {
		var index = 0;
		for(var i = 0; i < this.width; i++) {
			for(var j = 0; j < this.height; j++) {
				index = FLUID_IX(i, j);
				this.r[index] = this.rOld[index] = Math.random();
				if(this.isRGB) {
					this.g[index] = this.gOld[index] = Math.random();
					this.b[index] = this.bOld[index] = Math.random();
				}
			} 
		}
	}	

	this.getIndexForCellPosition = function(i, j) {
		if(i < 1) i=1; else if(i > _NX) i = _NX;
		if(j < 1) j=1; else if(j > _NY) j = _NY;
		return FLUID_IX(i, j);
	}

	this.getIndexForNormalizedPosition = function(x, y) {
		return this.getIndexForCellPosition(parseInt(x * _NX2), parseInt(y * _NY2));
	}

}