function particle(x, y, mass, window, solver) {
	
	this.x = x;
	this.y = y;
	this.mass = mass;
	this.vx = 0.0;
	this.vy = 0.0;
	this.width = window.width;
	this.height = window.height;
	this.solver = solver;
	
	this.addForce = function(fx, fy) {
//		console.log(fx + " < " + this.mass);
		if (Math.abs(fx)>this.mass) this.vx += fx/this.mass;
		if (Math.abs(fy)>this.mass) this.vy += fy/this.mass;		
	}
	
	this.update = function() {
		
		//use matrix to avg out surrounding u/vs

		
		
//		console.log(this.vx + ", " + this.vy);
		this.x += this.vx;
		this.y += this.vy;
		if (x<0) x = 0;
		if (x>this.width) x = this.width;
		if (y<0) y = 0;
		if (y>this.height) y = this.height;
		this.vx/=3.0;
		this.vy/=3.0;
	}	
}