var Engine = Matter.Engine;
var World = Matter.World;
var Bodies = Matter.Bodies;
var Composites = Matter.Composites;
var Composite = Matter.Composite;
var Constraint = Matter.Constraint;
var Sleeping = Matter.Sleeping;
var engine;
var bodies;
var world;

var molID;
var molarray = [];
var gcoll = 0;
var totalRows;

//store NaCl Properties
var Na;
var Cl;
var NaParticles = [];
var ClParticles = [];
var constraints = [];

//Store CaCl Properties
var Ca;
var Chl;
var Chl2;
var CaParticles = [];
var ChlParticles = [];
var Chl2Particles = [];
var constraints = [];

//Simulation Titles
var pageTitle=document.title;
var sim1="Introduction to Solutions";
var sim2="Solubility";
var sim3="Connected Chemistry - Physical Changes";
var sim4="Connected Chemistry - Chemical changes";
var sim5="";

//store propeties of bodies
var bodynumber = [];
var bodyPositionX = [];
var bodyPositionY = [];
var bodyType = [];
var minTemp = -20;
var bodycolgrnd = [];
var Ke = 0;
var TotalKE = 0;
var test = [];
var img;
var systemtemp;
var val=0;
var randomnum;
var mol1;
var mol2;
var idealRoomTemp = 25.00;
var currRoomTemp = idealRoomTemp;
var intVal;
var OneNegUnit = 8.8;
var OnePosUnit = 47.70;
var RateOfChangeNeg = OneNegUnit / 40;
var RateOfChangePos = OnePosUnit / 40;
var currRateChange = 0;
var collisionHeight = 0;
var totalSeconds = 0;
var timerVar;
var totalMolecules=0;
var totalMolecules1=0;
var totalMolecules2=0;
var stack = [];
var stackCounter = 0;
var radius=20;
var radius1=20;
var radius2=20;
var moleculenum=28;
var row,column;
var compounds = database['compounds'];
var freezingPoint;
var boilingPoint;




// loading IMAGE for different options
function setMolecule(ID) {

    stackCounter=0;
    if (ID == 1) {
        img = loadImage('Water.png');
        img1 = loadImage('Sodium-I.png');
        img2 = loadImage('Chloride.png');
        document.getElementById("molimg").src="Sodium-Chloride.png";
        document.getElementById("molname").innerHTML="NaCl";
        window.myLine.data.datasets[1].label="Sodium Chloride";
        setFreezingBoilingPoint("Water");
        radius1= 15;
        radius2 = 20;
        moleculenum=28;
        molID=1;
    }  else if (ID == 3) {
        setup();
        img = loadImage('Water.png');
        img1 = loadImage('Calcium-Ion.png');
        img2 = loadImage('Chloride.png');
        document.getElementById("molimg").src="Calcium-Chloride.png";
        document.getElementById("molname").innerHTML="Calcium-Chloride";
        window.myLine.data.datasets[1].label="Calcium-Chloride";
        setFreezingBoilingPoint("Water");
        moleculenum=28;
        molID=3;
    } else if (ID == 6) {
        setup();
        img = loadImage('Water.png');
        img1 = loadImage('Bicarbonate.png');
        img2 = loadImage('Sodium-I.png');
        document.getElementById("molimg").src="Sodium-Bicarbonate.png";
        document.getElementById("molname").innerHTML="Sodium Bicarbonate";
        window.myLine.data.datasets[1].label="Sodium Bicarbonate";
        setFreezingBoilingPoint("Water");
        radius1= 25;
        radius2 = 15;
        moleculenum=28;
        molID=6;
    } 
    setup();
} 

function setFreezingBoilingPoint(name){
    
    for(i = 0 ; i < compounds.length ; i++) //array to check each record
    {
        if(compounds[i].name == name)
        {
            freezingPoint=compounds[i].freezingPointCelsius;
            boilingPoint=compounds[i].boilingPointCelsius;
        }
    }
}

function preload(){

    img = loadImage('Water.png');
    molID = 1;
    setFreezingBoilingPoint("Water");
    if (pageTitle==sim1) {
        img = loadImage('Water.png');
        img1 = loadImage('Sodium-I.png');
        img2 = loadImage('Chloride.png');
        sliderHeat.disable();
    }
    else if(pageTitle==sim2){
        img = loadImage('Water.png');
        img1 = loadImage('Sodium-I.png');
        img2 = loadImage('Chloride.png');
        addDataset("NaCl");
    }
}

function setup() {

    imageMode(CENTER);
    frameRate(15);

    //shasaboo
    stackCounter=0;

    //Create and position the canvas
    canvasHeight = window.innerHeight*0.85;
    canvasWidth = window.innerWidth*0.7;
    
    cnv = createCanvas(canvasWidth, canvasHeight);
    cnv.position(0.3*window.innerWidth, 0.15*window.innerHeight);
    
    //Setup water molecules by default
    molID == 1
    img = loadImage('Water.png');
    
    //create an engine
    engine = Engine.create();
    world = engine.world;

    //set gravity parameters
    var gravityX = world.gravity.x;
    var gravityY = world.gravity.y;

    //add rectangles/floors
    var params1 = {
        isStatic: true,
        restitution: 1.0
    }

    var params = {
        isStatic: true
    }
    var ground = Bodies.rectangle(width / 2, height, width, 0.04*height*2, params1);
    var wall1 = Bodies.rectangle(0, height / 2, 0.016*width*2, height, params);
    var wall2 = Bodies.rectangle(width, height / 2, 0.016*width*2, height, params);
    var top = Bodies.rectangle(width / 2, 0, width, 0.025*height*2, params);
    World.add(world, [ground, wall1, wall2, top]);
    
       //reset bodies to 
    
    //add circle stack
   // addMolecules(28); 

    //RESET
    //slider
    document.getElementById("heatValue").innerHTML=0;
    sliderHeat.setValue(0);

    //timer
    totalSeconds=0;
    clearInterval(timerVar);
    var hour = Math.floor(totalSeconds /3600);
    var minute = Math.floor((totalSeconds - hour*3600)/60);
    seconds = totalSeconds - (hour*3600 + minute*60);
    timervalue = hour + ":" + minute + ":" + seconds;
    document.getElementById("timer").innerHTML = timervalue;

    //Graph
    removeData();

    //Room Temperature
    currRoomTemp = idealRoomTemp;
    document.getElementById("tempval").innerHTML = currRoomTemp.toFixed(3);

    //reset bodies to 0
    bodies = [];
    bodynumber = [];
    bodyType = [];
    NaParticles = [];
    ClParticles = [];
    console.log(bodyType);
    
    CaParticles = [];
    ChlParticles = [];
    Chl2Particles = [];

    //reset total molecules
    totalMolecules=0;
    totalMolecules1 = 0;
    //add circle stack
    addMolecules(moleculenum,molID);
    addMoleculesCacl(4,molID);
    //console.log(NaParticles);
    
}

function draw() {
    
    background(70);
    var timescale=document.getElementById("speedValue").innerHTML;
    drawinitialsetup();
    setHeatLine();
    Engine.update(engine);

    //heat slider value
    val = document.getElementById("heatValue").innerHTML;
    var intValtemp = currRoomTemp;

    //heatslider changes with value
    if (intValtemp > boilingPoint) {
        world.gravity.y = 0;
    } else if (intValtemp < freezingPoint) {
        world.gravity.y = 0.3;
        applyforce(0.005);
    } else if (intValtemp < boilingPoint) {
        world.gravity.y = 0.3;
        applyforce(0.00005);
    }

   applyforceNaCl(0.005);
    
    //Applying attrcative forces between all molecules
    function applyforce(mag) {
        for (var i = 0; i < bodies.length; i++) {
            for (var j = 0; j < bodies.length; j++) {
                var pvector = createVector((bodynumber[i].position.x), (bodynumber[i].position.y));
                var p2vector = createVector((bodynumber[j].position.x), (bodynumber[j].position.y));
                var distforce = dist(bodynumber[i].position.x, bodynumber[i].position.y, bodynumber[j].position.x, bodynumber[j].position.y);
                var diffvector = p5.Vector.sub(pvector, p2vector);
                diffvector.normalize();
                diffvector.mult(mag);
                //line(diffvector);
                if (distforce > 40) {
                    Matter.Body.applyForce(bodynumber[j], p2vector, diffvector);
                }
            }
        }
    }
      
    //Applying attractive forces between NaCl molecules
   function applyforceNaCl(mag) {
        for (var i = 0; i < NaParticles.length; i++) {
            for (var j = 0; j < NaParticles.length; j++) {
                var pvector = createVector((NaParticles[i].position.x),(NaParticles[i].position.y));
                var p2vector = createVector((NaParticles[j].position.x),(NaParticles[j].position.y));
                var distforce = dist(NaParticles[i].position.x, NaParticles[i].position.y, NaParticles[j].position.x, NaParticles[j].position.y);
                var diffvector = p5.Vector.sub(pvector, p2vector);
                diffvector.normalize();
                diffvector.mult(mag);
                if (distforce > 40) {
                    Matter.Body.applyForce(NaParticles[j], p2vector, diffvector);
                }
            }
        }
        for (var i = 0; i < ClParticles.length; i++) {
            for (var j = 0; j < ClParticles.length; j++) {
                var pvector = createVector((ClParticles[i].position.x),(ClParticles[i].position.y));
                var p2vector = createVector((ClParticles[j].position.x),(ClParticles[j].position.y));
                var distforce = dist(ClParticles[i].position.x, ClParticles[i].position.y, ClParticles[j].position.x, ClParticles[j].position.y);
                var diffvector = p5.Vector.sub(pvector, p2vector);
                diffvector.normalize();
                diffvector.mult(mag);
                if (distforce > 40) {
                    Matter.Body.applyForce(ClParticles[j], p2vector, diffvector);
                }
            }
        }
    }
    
    
    
   
    
      
    for (i = 0; i < bodies.length; i++) {
        collides(bodynumber[i]);
    }

    function collides(circlebody) { //implements ground collision

        collisionHeight = (circlebody.position.y + 20);

        if (collisionHeight >= (height - 35) && circlebody.collideFlag == -1) {

            gcoll += 1;
            circlebody.collideFlag = 1;

            intVal = document.getElementById("heatValue").innerHTML;

            if (intVal == 0) {
                intVal = idealRoomTemp;
            } else {
                intVal = intVal < 0 ? intVal * OneNegUnit + idealRoomTemp : intVal * OnePosUnit + idealRoomTemp;
            }

            //set current rate change according to currRoomTemp
            currRateChange = currRoomTemp < idealRoomTemp ? RateOfChangeNeg : RateOfChangePos;

            if (currRoomTemp > intVal) {
                currRoomTemp = currRoomTemp - currRateChange;
                if (currRoomTemp < intVal) {
                    currRoomTemp = intVal;
                }
            }
            if (currRoomTemp < intVal) {
                currRoomTemp = currRoomTemp + currRateChange;
                if (currRoomTemp > intVal) {
                    currRoomTemp = intVal;
                }
            }

            //Add temperature value onscreen
            document.getElementById("tempval").innerHTML = currRoomTemp.toFixed(3);

            bodycolgrnd.push(circlebody);
            if (val == 0) {
                circlebody.velocity.y = -abs(circlebody.velocity.x);
            } else if (val > 0) {
                circlebody.velocity.y = -2 * abs(circlebody.velocity.x);
                Ke += Ke;
            } else if (val < 0) {
                Ke -= Ke;
            }
        } else {

            if ((circlebody.position.y + 20) < (height - 35)) {
                circlebody.collideFlag = -1;
            } else // if collision not detected only due to the reason it was already in collided positon 
            {
                if (intVal < idealRoomTemp) // special case to deal for negative slider position, because molecules will remain in touch with bar
                {
                    //Probability for detection should be managed in case of negative slider position
                    if ((Math.random() * (11 - 1) + 1) < 3) // generate random num b/w 1 and 10
                    {
                        circlebody.collideFlag = -1;
                    }
                }
            }
        }
    }

    //set forces after and before separation
    if (gcoll == 0) {
        applyforce(0.0008);
        var mag = 0.0008;
    } else {
        applyforce(0);
        mag = 0;
    }

    //check circle collision
    function intersects(first, other) {
        var test=false;
        var d = dist(first.position.x, first.position.y, other.position.x, other.position.y);
        if (d <= 56) {
            first.velocity.x = -abs(first.velocity.x);
            first.velocity.y = -abs(first.velocity.y);
            other.velocity.x = -abs(other.velocity.x);
            other.velocity.y = -abs(other.velocity.y);
            test = true;
        }
        return test;
    }

    //Apply force between molecules
    for (var i = 0; i < bodies.length; i++) {
        var positionvec = createVector((bodynumber[i].position.x), (bodynumber[i].position.y));
        var v = createVector((0.02 * random(-1, 1)), (0.05 * random(-1, 1)));
        Matter.Body.applyForce(bodynumber[i], positionvec, v);
    }
}


function drawinitialsetup() {
    
    drawcircle();

    //draw floor
    stroke(108, 108, 108);
    strokeWeight(10);
    
    //bottom wall
    fill(108, 108, 108);
    rect(0, height - 0.04*height, width, 0.04*height);

    //left wall
     fill(108, 108, 108);
     rect(0, 0, 0.016*width, height);

     //right wall
     fill(108, 108, 108);
     rect(width - 15, 0, 0.016*width, height);

     //top wall
     fill(108, 108, 108);
     rect(0, 0, width, 0.025*height);
}

function setHeatLine(){

    var heat=document.getElementById("heatValue").innerHTML;
    if(heat>0)
        stroke(255*(heat/5), 0, 0);
    else if(heat<0)
        stroke(0, 0, 255*(-heat/5));
    strokeWeight(10);
    line(25, height - 0.04*height, width-25, height - 0.04*height);
}

function drawcircle() {
    stroke(255);
    strokeWeight(2);
    fill(255, 50);
    for (var i = 0; i < bodies.length; i++) {

        var circle = bodies[i];
        var pos = circle.position;
        var r = 20;
        var angle = circle.angle;
        push();
        translate(pos.x, pos.y);
        rotate(angle);
        //if (molarray.indexOf(i) != -1)
        //{
            //image(img1, 0, 0, 40, 28);
        //} else {
        /*if(bodyType[i]==10)
            image(img1, 0, 0);
        else if(bodyType[i]==11)
            image(img2, 0, 0);
        else */
            image(img, 0, 0);
        //}
        pop();
    }
    
    //Drawing Calcium Chloride
    if(molID == 3)
    {
    stroke(255);
    strokeWeight(2);
    fill(255, 50);
   for (var i = 0; i < CaParticles.length;i++)
   {
    var Ca = CaParticles[i];
    var Chl = ChlParticles[i];
    var Chl2 = Chl2Particles[i];
    var Capos = Na.position;
    var Chlpos = Chl.position;
    var Chl2pos = Chl2.position;
    var Caangle = Na.angle;
    var Chlangle = Cl.angle;
    var Chl2angle = Cl2.angle;
        
        push();
        translate(Capos.x, Napos.y);
        //rotate(Naangle);
        image(img2, 0, 0,30,30);
        pop();
        
        push();
        translate(Chlpos.x, Chlpos.y);
        //rotate(Clangle);
        image(img1, 0, 0,40,40);
        pop();
        
         push();
        translate(Chl2pos.x, Chl2pos.y);
        //rotate(Naangle);
        image(img2, 0, 0,30,30);
        pop();
    } 
    }   
    
}
function setHeatSliderVal(){
    document.getElementById("heatValue").innerHTML=document.getElementById("ex4").value;
}
function setAddMoleculeLabelNacl(){
    document.getElementById("moleculevalNacl").innerHTML=document.getElementById("ex11").value;
}
function setSpeedSliderVal(){
    document.getElementById("speedValue").innerHTML=document.getElementById("ex2").value;
}
function setAddMoleculeLabel(){
    document.getElementById("moleculeval").innerHTML=document.getElementById("ex1").value;
}
function addMoleculesBtn(){
    addMolecules(document.getElementById("moleculeval").innerHTML,molID);
}
function addMoleculesBtnNacl(){
    addMoleculesCacl(document.getElementById("moleculevalNacl").innerHTML,molID);
}

function addMoleculesCacl(number,ID)
{
     function particleCa(x, y) {
        var params = {
            restitution: 1.0,
            friction: 0,
            offset: 0,
            mass: 25,
            collideFlag: -1

        }
        return Bodies.circle(x, y, radius1, params);
    }
    
    function particleChl(x, y) {
        var params = {
            restitution: 1.0,
            friction: 0,
            offset: 0,
            mass: 16,
            collideFlag: -1

        }
        return Bodies.circle(x, y, radius2, params);
    }
    function particleChl2(x, y) {
        var params = {
            restitution: 1.0,
            friction: 0,
            offset: 0,
            mass: 16,
            collideFlag: -1

        }
        return Bodies.circle(x, y, radius2, params);
    }

 //totalMolecules1 += parseInt(number);
  

   if(molID == 3)
   {
   var x=400,y=250;
   for (var i= 0; i< number;i++,x+=85)
   {
   CaParticles.push(new particleCa(x,y));
   ChlParticles.push(new particleChl(x+40,y));
   Chl2Particles.push(new particleChl2(x+80,y));
   

   
   if (CaParticles.length < 5)
    {
    CaParticles[i].isSleeping = true;
    ChlParticles[i].isSleeping = true;
    Chl2Particles[i].isSleeping = true;
    console.log("working");
    //console.log(NaParticles[i].isSleeping);
    }
   }
   }
   
    //console.log(NaParticles.length + "  NA  " +ClParticles.length + "Cl");
    
    for (var i = (CaParticles.length - number)  ; i < CaParticles.length ;i++)
    {
    World.add(world, [CaParticles[i],ChlParticles[i],Chl2Particles[i]]);
    
    var options = {
    	bodyA : CaParticles[i],
    	bodyB : ChlParticles[i],
    	length : 40,
        stiffness : 1
        }
    var options1 = {
    	bodyA : ChlParticles[i],
    	bodyB : Chl2Particles[i],
    	length : 40,
        stiffness : 1
        }
    

    var constraint1 = Constraint.create(options);
    var constraint2 = Constraint.create(options1);
    
    World.add(world,[constraint1,constraint2]);
    constraints[i] = constraint;
    console.log(constraints[i]);
    }
}


////////////////for CaCl
var particleOptions = { 
        friction: 0.05,
        frictionStatic: 0.1,
        render: { visible: true } 
    };

var constraintOptions = { 
        render: { visible: false } 
    };

var softBody = Composites.softBody(450, 200, 10, 5, 0, 0, true, 15, particleOptions, constraintOptions);

World.add(world, [softBody]);

///////////////caCl


function setrowcolumn(number){

    if(number==28){
        row=7;
        column=4;
    }
    if(number==20){
        row=5;
        column=4;
    }
    if(number==12){
        row=4;
        column=3;
    }
    if(number==15){
        row=5;
        column=3;
    }
    if(number==40){
        row=8;
        column=5;
    }
    if (number<=10) {
        row=number/2;
        column=number/row;
    }

}

function addMolecules(number,ID){

    setrowcolumn(number);

    function makeCircle(x, y) {
        var params = {
            restitution: 1.0,
            friction: 0,
            offset: 0,
            mass: 18,
            collideFlag: -1
        }
        return Bodies.circle(x, y, radius, params);
    }

    // x, y, columns, rows, column gap, row gap
    stack[stackCounter] = Composites.stack(360, 20, row, column, 0, 0, makeCircle);
    totalMolecules+=row*column;
    World.add(world, stack[stackCounter]);
    var bodynum=bodynumber.length;
    if(bodies != undefined){
        var oldbodies=bodies;
        var newbodies= stack[stackCounter].bodies;
        bodies=oldbodies.concat(newbodies);
    }
    else{
        bodies = stack[stackCounter].bodies;
    }

    for (var i = 0; i < bodies.length; i++) {
        bodynumber[i] = bodies[i];
        if (stackCounter==0) {
            bodynumber[i].isSleeping = true;
        }
        bodyPositionX[i] = bodies[i].position.x;
        bodyPositionY[i] = bodies[i].position.y;
        if(bodyType[i]==undefined)
            bodyType[i] = ID;
    }
    stackCounter++;
}
//start simulation
function start(){
    for (var i = 0; i < bodynumber.length ; i++){
        bodynumber[i].isSleeping = false;
    }
       for (var i = 0; i < NaParticles.length ; i++)
    {
        NaParticles[i].isSleeping = false;
        ClParticles[i].isSleeping = false;
    /*if (NaParticles.length < 3)
    {
    NaParticles[i].isSleeping = false;
    ClParticles[i].isSleeping = false;
    }*/
    }
    timerVar = setInterval(countTimer, 1000);
}
//pause the simulation
function pause(){
    for (var i = 0; i < bodynumber.length ; i++){
        bodynumber[i].isSleeping = true;
    }
    clearInterval(timerVar);
    for (var i = 0; i < NaParticles.length ; i++)
    {
    NaParticles[i].isSleeping = true;
    ClParticles[i].isSleeping = true;
    console.log("working..");
    console.log(NaParticles[i].isSleeping);
    console.log(ClParticles[i].isSleeping);
    }
}

//Timer
function countTimer() {
    ++totalSeconds;
    var hour = Math.floor(totalSeconds /3600);
    var minute = Math.floor((totalSeconds - hour*3600)/60);
    seconds = totalSeconds - (hour*3600 + minute*60);
    timervalue = hour + ":" + minute + ":" + seconds;
    document.getElementById("timer").innerHTML = timervalue;
    //document.getElementById("timer").style.color = "black";
    document.getElementById("timer").style.font = " bold 20px arial,serif";

    //Setting the graph
    addLabel(timervalue);
    addData(0, totalMolecules);
    addData(1, totalMolecules1);
}
