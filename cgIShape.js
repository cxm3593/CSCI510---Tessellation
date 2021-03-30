//
// fill in code that creates the triangles for a cube with dimensions 1x1x1
// on each side (and the origin in the center of the cube). with an equal
// number of subdivisions along each cube face as given by the parameter
//subdivisions
//
function makeCube (subdivisions)  {
    
    // fill in your code here.
    // delete the code below first.

	//// Old code
    //addTriangle (-0.5, 0.5, 0.5, 0.5, -0.5, 0.5, 0.5, 0.5, 0.5);

	// var A = [-0.5,  0.5,  0.5];
	// var B = [ 0.5,  0.5,  0.5];
	// var C = [-0.5, -0.5,  0.5];
	// var D = [ 0.5, -0.5,  0.5];
	// var E = [-0.5,  0.5, -0.5];
	// var F = [ 0.5,  0.5, -0.5];
	// var G = [-0.5, -0.5, -0.5];
	// var H = [ 0.5, -0.5, -0.5];

	// var triangles = [
	// 	[A, C, B],
	// 	[B, C, D],
	// 	[E, A, F],
	// 	[F, A, B],
	// 	[E, G, A],
	// 	[A, G, C],
	// 	[B, D, F],
	// 	[F, D, H],
	// 	[F, H, E],
	// 	[E, H, G],
	// 	[H, D, G],
	// 	[G, D, C]

	// ];

	// // while subdivision is not finished , generate more triangles from last triangles array
	// while (subdivisions > 0){
	// 	tri_num = triangles.length;
	// 	for(var i=0; i<tri_num; i++){
	// 		var this_tri = triangles.shift();
	// 		var m1 = [
	// 			(this_tri[0][0] + this_tri[1][0]) / 2,
	// 			(this_tri[0][1] + this_tri[1][1]) / 2,
	// 			(this_tri[0][2] + this_tri[1][2]) / 2
	// 		];
	// 		var m2 = [
	// 			(this_tri[0][0] + this_tri[2][0]) / 2,
	// 			(this_tri[0][1] + this_tri[2][1]) / 2,
	// 			(this_tri[0][2] + this_tri[2][2]) / 2
	// 		];
	// 		var m3 = [
	// 			(this_tri[1][0] + this_tri[2][0]) / 2,
	// 			(this_tri[1][1] + this_tri[2][1]) / 2,
	// 			(this_tri[1][2] + this_tri[2][2]) / 2
	// 		];
	// 		var v0 = this_tri[0];
	// 		var v1 = this_tri[1];
	// 		var v2 = this_tri[2];

	// 		triangles.push([v0, m1, m2]);
	// 		triangles.push([m2, m1, m3]);
	// 		triangles.push([m1, v1, m3]);
	// 		triangles.push([m2, m3, v2]);
	// 	}

	// 	subdivisions --;

	// }
	//// Old code

	var triangles = [];
	const step = 1 / subdivisions;

	// Front
	for(var i=0; i<subdivisions; i++){
		for(var j=0;j<subdivisions;j++){
			var v0 = [-0.5 + i*step, -0.5 + j*step, 0.5];
			var v1 = [-0.5 + (i+1)*step, -0.5 + j*step, 0.5];
			var v2 = [-0.5 + i*step, -0.5 + (j+1)*step, 0.5];
			var v3 = [-0.5 + (i+1)*step, -0.5 + (j+1)*step, 0.5];

			triangles.push([v0, v3, v2]);
			triangles.push([v0, v1, v3]);
		}
	}

	// Left
	for(var i=0; i<subdivisions; i++){
		for(var j=0;j<subdivisions;j++){
			var v0 = [-0.5, -0.5+i*step, -0.5+j*step];
			var v1 = [-0.5, -0.5+(i+1)*step, -0.5+j*step];
			var v2 = [-0.5, -0.5+i*step, -0.5+(j+1)*step];
			var v3 = [-0.5, -0.5+(i+1)*step, -0.5+(j+1)*step];

			triangles.push([v0, v2, v3]);
			triangles.push([v0, v3, v1]);
		}
	}

	// Right
	for(var i=0; i<subdivisions; i++){
		for(var j=0;j<subdivisions;j++){
			var v0 = [0.5, -0.5+i*step, -0.5+j*step];
			var v1 = [0.5, -0.5+(i+1)*step, -0.5+j*step];
			var v2 = [0.5, -0.5+i*step, -0.5+(j+1)*step];
			var v3 = [0.5, -0.5+(i+1)*step, -0.5+(j+1)*step];

			triangles.push([v0, v3, v2]);
			triangles.push([v0, v1, v3]);
		}
	}

	// Top
	for(var i=0; i<subdivisions; i++){
		for(var j=0;j<subdivisions;j++){
			var v0 = [-0.5+i*step, 0.5, -0.5+j*step];
			var v1 = [-0.5+(i+1)*step, 0.5, -0.5+j*step];
			var v2 = [-0.5+i*step, 0.5, -0.5+(j+1)*step];
			var v3 = [-0.5+(i+1)*step, 0.5, -0.5+(j+1)*step];

			triangles.push([v0, v2, v3]);
			triangles.push([v0, v3, v1]);
		}
	}

	// Back
	for(var i=0; i<subdivisions; i++){
		for(var j=0;j<subdivisions;j++){
			var v0 = [-0.5 + i*step, -0.5 + j*step, -0.5];
			var v1 = [-0.5 + (i+1)*step, -0.5 + j*step, -0.5];
			var v2 = [-0.5 + i*step, -0.5 + (j+1)*step, -0.5];
			var v3 = [-0.5 + (i+1)*step, -0.5 + (j+1)*step, -0.5];

			triangles.push([v0, v2, v3]);
			triangles.push([v0, v3, v1]);
		}
	}

	// Buttom
	for(var i=0; i<subdivisions; i++){
		for(var j=0;j<subdivisions;j++){
			var v0 = [-0.5+i*step, -0.5, -0.5+j*step];
			var v1 = [-0.5+(i+1)*step, -0.5, -0.5+j*step];
			var v2 = [-0.5+i*step, -0.5, -0.5+(j+1)*step];
			var v3 = [-0.5+(i+1)*step, -0.5, -0.5+(j+1)*step];

			triangles.push([v0, v3, v2]);
			triangles.push([v0, v1, v3]);
		}
	}

	// add triangles to finish the make process
	for (tri of triangles){
		addTriangle(
			tri[0][0], tri[0][1], tri[0][2],
			tri[1][0], tri[1][1], tri[1][2],
			tri[2][0], tri[2][1], tri[2][2]
			);
	}

    
}


//
// fill in code that creates the triangles for a cylinder with diameter 1
// and height of 1 (centered at the origin) with the number of subdivisions
// around the base and top of the cylinder (given by radialdivision) and
// the number of subdivisions along the surface of the cylinder given by
//heightdivision.
//
function makeCylinder (radialdivision,heightdivision){
    // fill in your code here.
	var triangles = [];
	var circ_triangles = []; // Contains triangles from the top and bottom
	var side_triangles = []; // Contains triagnles from the side and 

	var Bottom_Center = [0, -0.5, 0];
	var Top_Center = [0, 0.5, 0];
	
	const radius = 0.5;

	
	var alpha_deg = 0.0;
	var step = 360 / radialdivision;
	var vertical_step = 1 / heightdivision;

	// Generate top and bottom triangles
	for(var i=0; i<radialdivision; i++){


		var b0 = [radius * Math.cos(radians(alpha_deg)), -0.5, radius * Math.sin(radians(alpha_deg))];
		var t0 = [radius * Math.cos(radians(alpha_deg)),  0.5, radius * Math.sin(radians(alpha_deg))];
		alpha_deg += step;
		var b1 = [radius * Math.cos(radians(alpha_deg)), -0.5, radius * Math.sin(radians(alpha_deg))];
		var t1 = [radius * Math.cos(radians(alpha_deg)),  0.5, radius * Math.sin(radians(alpha_deg))];

		circ_triangles.push([Bottom_Center, b0, b1]);
		circ_triangles.push([Top_Center, t1, t0]);

		// side_triangles.push([t0, b0, t1]);
		// side_triangles.push([t1, b0, b1]);

		var height_mark = -0.5;

		for(var j=0; j<heightdivision; j++){
			var m0 = [b0[0], height_mark, b0[2]];
			var m1 = [b1[0], height_mark, b1[2]];
			height_mark += vertical_step;
			var m2 = [b0[0], height_mark, b0[2]];
			var m3 = [b1[0], height_mark, b1[2]];

			side_triangles.push([m3, m1, m0]);
			side_triangles.push([m2, m3, m0]);
		}
	}

	// For Debug
	// console.log(radialdivision, " + ", circ_triangles.length);
	// for(d of circ_triangles){
	// 	console.log(d);
	// }

	triangles = circ_triangles.concat(side_triangles);

	// add triangles to finish the make process
	for (tri of triangles){
		addTriangle(
			tri[0][0], tri[0][1], tri[0][2],
			tri[1][0], tri[1][1], tri[1][2],
			tri[2][0], tri[2][1], tri[2][2]
			);
	}
	
}


//
// fill in code that creates the triangles for a cone with diameter 1
// and height of 1 (centered at the origin) with the number of
// subdivisions around the base of the cone (given by radialdivision)
// and the number of subdivisions along the surface of the cone
//given by heightdivision.
//
function makeCone (radialdivision, heightdivision) {
    // fill in your code here.
	var triangles = [];
	var circ_triangles = []; // Contains triangles from the top and bottom
	var side_triangles = []; // Contains triangles from the side and 

	var Bottom_Center = [0, -0.5, 0];
	var Apex = [0, 0.5, 0];

	const radius = 0.5;
	var alpha_deg = 0.0;
	var step = 360 / radialdivision;

	

	for(var i=0; i<radialdivision; i++){
		var b0 = [radius * Math.cos(radians(alpha_deg)), -0.5, radius * Math.sin(radians(alpha_deg))];
		alpha_deg += step;
		var b1 = [radius * Math.cos(radians(alpha_deg)), -0.5, radius * Math.sin(radians(alpha_deg))];

		circ_triangles.push([Bottom_Center, b0, b1]);

		//side_triangles.push([Apex, b1, b0]);
		var m0 = [b0[0], b0[1], b0[2]];
		var m1 = [b1[0], b1[1], b1[2]];
		var x0_step = (Apex[0] - b0[0]) / heightdivision;
		var y0_step = (Apex[1] - b0[1]) / heightdivision;
		var z0_step = (Apex[2] - b0[2]) / heightdivision;
		var x1_step = (Apex[0] - b1[0]) / heightdivision;
		var y1_step = (Apex[1] - b1[1]) / heightdivision;
		var z1_step = (Apex[2] - b1[2]) / heightdivision;

		for(var j=0; j<heightdivision; j++){
			if(j == heightdivision - 1){
				side_triangles.push([Apex, m1, m0]);
			}
			else{
				var m2 = [m0[0] + x0_step, m0[1]+y0_step, m0[2]+z0_step];
				var m3 = [m1[0] + x1_step, m1[1]+y1_step, m1[2]+z1_step];
				side_triangles.push([m2, m1, m0]);
				side_triangles.push([m3, m1, m2]);

				m0 = m2;
				m1 = m3;
			}
		}
	}

	triangles = circ_triangles.concat(side_triangles);

	// add triangles to finish the make process
	for (tri of triangles){
		addTriangle(
			tri[0][0], tri[0][1], tri[0][2],
			tri[1][0], tri[1][1], tri[1][2],
			tri[2][0], tri[2][1], tri[2][2]
			);
	}
}
    
//
// fill in code that creates the triangles for a sphere with diameter 1
// (centered at the origin) with number of slides (longitude) given by
// slices and the number of stacks (lattitude) given by stacks.
// For this function, you will implement the tessellation method based
// on spherical coordinates as described in the video (as opposed to the
//recursive subdivision method).
//
function makeSphere (slices, stacks) {
    // fill in your code here.
	var triangles = [];
	var origin = [0.0, 0.0, 0.0];
	const radius = 0.5;
	const longi_step = radians(360 / slices); // in radian
	const lati_step = radians(180 / stacks); // in radian

	var theta = 0.0;
	

	for (var i=0; i<slices; i++){

		var phi = 0.0;
		for(var j=0; j<stacks; j++){
			var v1 = [radius * Math.sin(theta) * Math.sin(phi), radius * Math.cos(phi), radius * Math.cos(theta) * Math.sin(phi)];
			var v2 = [radius * Math.sin(theta + longi_step) * Math.sin(phi), radius * Math.cos(phi), radius * Math.cos(theta + longi_step) * Math.sin(phi)];
			var v3 = [radius * Math.sin(theta) * Math.sin(phi + lati_step), radius * Math.cos(phi + lati_step), radius * Math.cos(theta) * Math.sin(phi + lati_step)];
			var v4 = [radius * Math.sin(theta + longi_step) * Math.sin(phi + lati_step), radius * Math.cos(phi + lati_step), radius * Math.cos(theta + longi_step) * Math.sin(phi + lati_step)];

			triangles.push([v1, v3, v2]);
			triangles.push([v2, v3, v4]);

			phi += lati_step;
		}

		theta += longi_step;
	}

	console.log(triangles.length);

	for (tri of triangles){
		addTriangle(
			tri[0][0], tri[0][1], tri[0][2],
			tri[1][0], tri[1][1], tri[1][2],
			tri[2][0], tri[2][1], tri[2][2]
			);
	}


}


////////////////////////////////////////////////////////////////////
//
//  Do not edit below this line
//
///////////////////////////////////////////////////////////////////

function radians(degrees)
{
  var pi = Math.PI;
  return degrees * (pi/180);
}

function addTriangle (x0,y0,z0,x1,y1,z1,x2,y2,z2) {

    
    var nverts = points.length / 4;
    
    // push first vertex
    points.push(x0);  bary.push (1.0);
    points.push(y0);  bary.push (0.0);
    points.push(z0);  bary.push (0.0);
    points.push(1.0);
    indices.push(nverts);
    nverts++;
    
    // push second vertex
    points.push(x1); bary.push (0.0);
    points.push(y1); bary.push (1.0);
    points.push(z1); bary.push (0.0);
    points.push(1.0);
    indices.push(nverts);
    nverts++
    
    // push third vertex
    points.push(x2); bary.push (0.0);
    points.push(y2); bary.push (0.0);
    points.push(z2); bary.push (1.0);
    points.push(1.0);
    indices.push(nverts);
    nverts++;
}

