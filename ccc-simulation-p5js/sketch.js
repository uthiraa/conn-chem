var moleculeArray = [];

var imageDir = "img/svg/";

var numberOfParticles = 20;

function preload() {
}

function setup() {
    // Creates drawing canvas for P5 objects
    createCanvas(800, 600);

    // Adjust framerate to slow animation, useful for debugging
    // frameRate(5);

    // Sets all x,y coordinates to the center of images
    imageMode(CENTER);

    // Sets the line thickness of shapes
    strokeWeight(3);

    // Calculates how many rows of particles to draw, currently assumes 5 columns of particles
    var rowsOfParticles = Math.floor(numberOfParticles / 5)

    // Increments the number of rows if there isn't exactly a multiple of 5 particles
    if (numberOfParticles % 5 > 0) {
        rowsOfParticles++;
    }

    // Creates molecule collection in a grid with random velocity vectors
    for (var i = 0; i < rowsOfParticles; i++) {
        for (var j = 0; j < 5; j++) {
            moleculeArray[5 * i + j] = new Particle(161, createVector(width / 10 + (j * (width / 5)), height / 10 + i * (height / 10)), p5.Vector.random2D().mult(random(1,5)), 5 * i + j);
        }
    }

    // Creates molecule collection with random velocity vectors
    // More advanced behavior is handled above, this only creates one row of particles
    // for (var i = 0; i < numberOfParticles; i++) {
    //     moleculeArray[i] = new Particle(25, createVector((i + 1) * width / 5 - width / 10, height / 4), p5.Vector.random2D().mult(4), i);
    // }

}

function draw() {
    // The background must be redrawn with each draw loop to avoid molecules leaving permanent traces
    background(204, 204, 204);

    // Kinetic energy reset after each loop
    var systemKineticEnergy = 0;

    // Array of particles - this is iterated on each draw loop to render the objects to the canvas
    for (var i = 0; i < numberOfParticles; i++) {
        moleculeArray[i].drawParticle();
        moleculeArray[i].checkCollision();

        for (var j = 0; j < numberOfParticles; j++) {
            if ((moleculeArray[i].id == moleculeArray[j].id) || moleculeArray[j].indexed) {
                // console.log("skipped");
                continue;
            } else {
                // console.log("checking for collision");
                moleculeArray[i].checkParticleCollision(moleculeArray[j]);
            }
        }

        moleculeArray[i].updatePosition();
        moleculeArray[i].updateKineticEnergy();

        // moleculeArray[i].showLabel();
        // moleculeArray[i].showVelocity();
        // moleculeArray[i].log();
    }

    // Loop to update properties and variables after positions and velocities have been calculated
    for (var i = 0; i < numberOfParticles; i++) {
        // Resets the boolean for whether a particular molecule was already iterated over
        // This is is meant to trim down the number of collision calculations being performed each loop
        moleculeArray[i].indexed = false;

        // Get KE of each particle
        systemKineticEnergy += moleculeArray[i].getKineticEnergy();
    }

    // Display KE of system
    // text(systemKineticEnergy, width / 2, height - 20);

}

// Particle class to handle each particle and its properties
class Particle {

    constructor(key, position, velocity, id) {
        // This is the value in the "ID" field of the database variable
        this.databaseKey = key;

        // This gives each molecule a unique index to track it within a set of molecules
        this.id = id;

        // Particle ID is database key, used to pull name
        this.name = database[this.databaseKey - 1].name;

        // Sets the position of the particle
        this.position = createVector(position.x, position.y);

        // Sets the velocity of the particle
        this.velocity = createVector(velocity.x, velocity.y)

        // Initializes the particle mass
        this.mass = (random() < 0.5) ? 1 : 2;

        // Sets the KE of the particle
        this.kineticEnergy = this.getKineticEnergy();

        // Sets the image reference for using the ID number of the particle
        this.imageUrl = imageDir + database[this.databaseKey - 1].file;
        this.imageObject = loadImage(this.imageUrl);

        // Collision is a class that allows you to reference which particles were in a collision, useful for post collision behavior
        this.collision = null;

        // Keeps track of whether a molecule was already indexed in a pass of velocity calculations. 
        // Velocity is already set by partner particle, so recalculating it with the same method would give the wrong result
        this.indexed = false;

    }

    /**
     * Collision detection functions
     */

    // Checks whether the particle collides with the boundary
    checkCollision() {
        // Calculates the width of the particle glyph
        var imageWidth = this.imageObject.width;
        var imageHeight = this.imageObject.height;

        // Checks if width of object is within horizontal collision distance of boundary
        if (this.position.x < imageWidth / 2.0 || (width - this.position.x) < imageWidth / 2.0) {

            // This logic handles the case where particles get stuck to the wall
            if (this.position.x > width / 2) {
                this.velocity.x = -abs(this.velocity.x);
            } else {
                this.velocity.x = abs(this.velocity.x);
            }
        }

        // Check if height of object is within vertical collision distance of boundary
        if (this.position.y < imageHeight / 2.0 || (height - this.position.y) < imageHeight / 2.0) {

            // This logic handles the case where particles get stuck to the wall
            if (this.position.y > height / 2) {
                this.velocity.y = -abs(this.velocity.y);
            } else {
                this.velocity.y = abs(this.velocity.y);
            }
        }
    }

    // Checks whether two particles are going to collide
    checkParticleCollision(particle) {
        var distanceBetweenParticles = p5.Vector.sub(particle.position, this.position).mag();
        var collisionDistance = this.getRadius() + particle.getRadius();

        // Evaluates whether the particles are within a collision radius of one another
        if (distanceBetweenParticles <= collisionDistance) {

            // If there is currently no collision object stored on this class, that means the particle could collide
            if (this.collision == null) {

                // Updating this particle
                this.collision = new Collision(this, particle);
                var dv1 = this.collision.calculateVelocityVersion2();

                // Updating other particle
                particle.collision = new Collision(particle, this);
                var dv2 = particle.collision.calculateVelocityVersion2();

                // Updating both velocities, note that we must find the value of the velocities first before updating since the calculations of both require the initial state values
                this.velocity.add(dv1);
                particle.velocity.add(dv2);

                particle.indexed = true;

            }

        }
        // This branch is very important because it prevents velocities from updating while particles are still within collision distance of one another
        // Without setting this, the particle behavior after collision was very erratic
        else if (distanceBetweenParticles > collisionDistance) {
            // This checks if the other particle is indexed in the current collision object
            if (this.collision != null && this.collision.contains(particle)) {
                this.collision = null;
                particle.collision = null;
            }
        }


    }

    /**
     * Calculate properties of particle
     */

    // Updates the x,y coordinates of the particle
    updatePosition() {
        this.position.add(this.velocity);
    }

    updateKineticEnergy() {
        this.kineticEnergy = 0.5 * this.mass * Math.pow(this.velocity.mag(), 2);
    }

    /**
     * Get functions
     */

    // Returns a P5 Vector containing the particle velocity
    getVelocity() {
        return this.velocity;
    }

    // Provides logic to determine the collision radius of the particle
    getRadius() {
        // Returns the max value of the image width, this will need to be tweaked for long molecules like pentane
        return (this.imageObject.width >= this.imageObject.height) ? this.imageObject.width / 2.0 : this.imageObject.height / 2.0;
    }

    getKineticEnergy() {
        return this.kineticEnergy;
    }

    /**
     * Set functions
     */

    // Updates the vx,vy based on the force field
    setVelocity(velocity) {
        this.velocity.x = velocity.x;
        this.velocity.y = velocity.y;
    }


    /*
    * Rendering functions
    */

    // Renders the image glyph to canvas using P5.js Image object
    drawParticle() {
        image(this.imageObject, this.position.x, this.position.y);
    }

    // Displays a velocity vector representation for each particle
    showVelocity() {

        line(this.position.x, this.position.y, (this.position.x + 3 * this.velocity.x), (this.position.y + 3 * this.velocity.y));
    }

    showLabel() {
        // Displays a hovering bookeeping index for each particle
        text(this.id, this.position.x + 25, this.position.y + 25);
    }

    /**
     * Logging functions
    */

    // Writes particle ID, position, and velocity to the console
    log() {
        console.log(this.id + "\t" + this.position.x + "\t" + this.position.y + "\t" + this.velocity.x + "\t" + this.velocity.y);
    }

    // Provides a string representation of the Particle object
    toString() {
        return "This particle has the chemical identity '" + this.name + "'";
    }

};

class Collision {
    constructor(particle1, particle2) {
        this.particle1 = particle1;
        this.particle2 = particle2;
    }

    // Calculates the post collision velocity 
    calculateVelocity() {

        var reducedMass = (2.0 * this.particle2.mass) / (this.particle1.mass + this.particle2.mass);
        var dr = p5.Vector.sub(this.particle2.position, this.particle1.position);
        dr.normalize();
        var dot = dr.dot(this.particle1.velocity) * -1.0 * reducedMass;
        dr.mult(dot);

        return dr;
    }

    calculateVelocityVersion2() {
        // var initialVelocity = createVector(this.velocity.x, this.velocity.y);
        var reducedMass = (2.0 * this.particle2.mass) / (this.particle1.mass + this.particle2.mass);
        var dx = p5.Vector.sub(this.particle1.position, this.particle2.position);
        var dxmag = dx.mag();
        var dv = p5.Vector.sub(this.particle1.velocity, this.particle2.velocity);

        var dotProduct = p5.Vector.dot(dv, dx);

        return dx.mult((-reducedMass * dotProduct) / Math.pow(dxmag, 2));
    }

    // Returns true if the particle in the argument is one of the particles involved in the collision
    contains(particle) {
        return (particle.id == this.particle1.id || particle.id == this.particle2.id) ? true : false;
    }

    toString() {
        return "This collision contains particle " + this.particle1.id + " and particle " + this.particle2.id;
    }
};

// Array of Javascript objects containing names and file locations for different species
var database = [
    {
        "ID": 1,
        "file": "384100F6.png",
        "name": "384100F6",
        "nameSlug": "384100f6",
        "fileExtension": "png"
    },
    {
        "ID": 2,
        "file": "6161D219.png",
        "name": "6161D219",
        "nameSlug": "6161d219",
        "fileExtension": "png"
    },
    {
        "ID": 3,
        "file": "6161D21A.png",
        "name": "6161D21A",
        "nameSlug": "6161d21a",
        "fileExtension": "png"
    },
    {
        "ID": 4,
        "file": "Acetate.svg",
        "name": "Acetate",
        "nameSlug": "acetate",
        "fileExtension": "svg"
    },
    {
        "ID": 5,
        "file": "Acetic-Acid.svg",
        "name": "Acetic Acid",
        "nameSlug": "acetic-acid",
        "fileExtension": "svg"
    },
    {
        "ID": 6,
        "file": "Aluminum-Chloride.svg",
        "name": "Aluminum Chloride",
        "nameSlug": "aluminum-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 7,
        "file": "Aluminum-Ion.svg",
        "name": "Aluminum Ion",
        "nameSlug": "aluminum-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 8,
        "file": "Aluminum-Oxide.svg",
        "name": "Aluminum Oxide",
        "nameSlug": "aluminum-oxide",
        "fileExtension": "svg"
    },
    {
        "ID": 9,
        "file": "Aluminum.svg",
        "name": "Aluminum",
        "nameSlug": "aluminum",
        "fileExtension": "svg"
    },
    {
        "ID": 10,
        "file": "Ammonia-Dots.svg",
        "name": "Ammonia Dots",
        "nameSlug": "ammonia-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 11,
        "file": "Ammonia.svg",
        "name": "Ammonia",
        "nameSlug": "ammonia",
        "fileExtension": "svg"
    },
    {
        "ID": 12,
        "file": "Ammonium-Acetate.svg",
        "name": "Ammonium Acetate",
        "nameSlug": "ammonium-acetate",
        "fileExtension": "svg"
    },
    {
        "ID": 13,
        "file": "Ammonium-Chloride.svg",
        "name": "Ammonium Chloride",
        "nameSlug": "ammonium-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 14,
        "file": "Ammonium-Dots.svg",
        "name": "Ammonium Dots",
        "nameSlug": "ammonium-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 15,
        "file": "Ammonium-Nitrate.svg",
        "name": "Ammonium Nitrate",
        "nameSlug": "ammonium-nitrate",
        "fileExtension": "svg"
    },
    {
        "ID": 16,
        "file": "Ammonium.svg",
        "name": "Ammonium",
        "nameSlug": "ammonium",
        "fileExtension": "svg"
    },
    {
        "ID": 17,
        "file": "Barium-144.svg",
        "name": "Barium 144",
        "nameSlug": "barium-144",
        "fileExtension": "svg"
    },
    {
        "ID": 18,
        "file": "Bicarbonate.svg",
        "name": "Bicarbonate",
        "nameSlug": "bicarbonate",
        "fileExtension": "svg"
    },
    {
        "ID": 19,
        "file": "Blue.svg",
        "name": "Blue",
        "nameSlug": "blue",
        "fileExtension": "svg"
    },
    {
        "ID": 20,
        "file": "Boron-Ion.svg",
        "name": "Boron Ion",
        "nameSlug": "boron-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 21,
        "file": "Boron.svg",
        "name": "Boron",
        "nameSlug": "boron",
        "fileExtension": "svg"
    },
    {
        "ID": 22,
        "file": "Boron-Tetrachloride-Dots.svg",
        "name": "Boron Tetrachloride Dots",
        "nameSlug": "boron-tetrachloride-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 23,
        "file": "Boron-Tetrachloride.svg",
        "name": "Boron Tetrachloride",
        "nameSlug": "boron-tetrachloride",
        "fileExtension": "svg"
    },
    {
        "ID": 24,
        "file": "Boron-Trichloride-Dots.svg",
        "name": "Boron Trichloride Dots",
        "nameSlug": "boron-trichloride-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 25,
        "file": "Boron-Trichloride.svg",
        "name": "Boron Trichloride",
        "nameSlug": "boron-trichloride",
        "fileExtension": "svg"
    },
    {
        "ID": 26,
        "file": "Bromide-Dots.svg",
        "name": "Bromide Dots",
        "nameSlug": "bromide-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 27,
        "file": "Bromide.svg",
        "name": "Bromide",
        "nameSlug": "bromide",
        "fileExtension": "svg"
    },
    {
        "ID": 28,
        "file": "Bromine-Atom.svg",
        "name": "Bromine Atom",
        "nameSlug": "bromine-atom",
        "fileExtension": "svg"
    },
    {
        "ID": 29,
        "file": "Bromine-Ion.svg",
        "name": "Bromine Ion",
        "nameSlug": "bromine-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 30,
        "file": "Bromine.svg",
        "name": "Bromine",
        "nameSlug": "bromine",
        "fileExtension": "svg"
    },
    {
        "ID": 31,
        "file": "Butane.svg",
        "name": "Butane",
        "nameSlug": "butane",
        "fileExtension": "svg"
    },
    {
        "ID": 32,
        "file": "Butene.svg",
        "name": "Butene",
        "nameSlug": "butene",
        "fileExtension": "svg"
    },
    {
        "ID": 33,
        "file": "Calcium-Chloride.svg",
        "name": "Calcium Chloride",
        "nameSlug": "calcium-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 34,
        "file": "Calcium-Ion.svg",
        "name": "Calcium Ion",
        "nameSlug": "calcium-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 35,
        "file": "Calcium.svg",
        "name": "Calcium",
        "nameSlug": "calcium",
        "fileExtension": "svg"
    },
    {
        "ID": 36,
        "file": "Carbide.svg",
        "name": "Carbide",
        "nameSlug": "carbide",
        "fileExtension": "svg"
    },
    {
        "ID": 37,
        "file": "Carbonate.svg",
        "name": "Carbonate",
        "nameSlug": "carbonate",
        "fileExtension": "svg"
    },
    {
        "ID": 38,
        "file": "Carbon-Dioxide.svg",
        "name": "Carbon Dioxide",
        "nameSlug": "carbon-dioxide",
        "fileExtension": "svg"
    },
    {
        "ID": 39,
        "file": "Carbonic-Acid.svg",
        "name": "Carbonic Acid",
        "nameSlug": "carbonic-acid",
        "fileExtension": "svg"
    },
    {
        "ID": 40,
        "file": "Carbon-Ion.svg",
        "name": "Carbon Ion",
        "nameSlug": "carbon-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 41,
        "file": "Carbon-Monoxide.svg",
        "name": "Carbon Monoxide",
        "nameSlug": "carbon-monoxide",
        "fileExtension": "svg"
    },
    {
        "ID": 42,
        "file": "Carbon.svg",
        "name": "Carbon",
        "nameSlug": "carbon",
        "fileExtension": "svg"
    },
    {
        "ID": 43,
        "file": "Catalyst.svg",
        "name": "Catalyst",
        "nameSlug": "catalyst",
        "fileExtension": "svg"
    },
    {
        "ID": 44,
        "file": "Cesium-137.svg",
        "name": "Cesium 137",
        "nameSlug": "cesium-137",
        "fileExtension": "svg"
    },
    {
        "ID": 45,
        "file": "Chloride-Dots.svg",
        "name": "Chloride Dots",
        "nameSlug": "chloride-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 46,
        "file": "Chloride.svg",
        "name": "Chloride",
        "nameSlug": "chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 47,
        "file": "Chlorine-45.svg",
        "name": "Chlorine 45",
        "nameSlug": "chlorine-45",
        "fileExtension": "svg"
    },
    {
        "ID": 48,
        "file": "Chlorine-Atom.svg",
        "name": "Chlorine Atom",
        "nameSlug": "chlorine-atom",
        "fileExtension": "svg"
    },
    {
        "ID": 49,
        "file": "Chlorine-Ion.svg",
        "name": "Chlorine Ion",
        "nameSlug": "chlorine-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 50,
        "file": "Chlorine-Nitrate.svg",
        "name": "Chlorine Nitrate",
        "nameSlug": "chlorine-nitrate",
        "fileExtension": "svg"
    },
    {
        "ID": 51,
        "file": "Chlorine.svg",
        "name": "Chlorine",
        "nameSlug": "chlorine",
        "fileExtension": "svg"
    },
    {
        "ID": 52,
        "file": "Chromium-III-Oxide.svg",
        "name": "Chromium III Oxide",
        "nameSlug": "chromium-iii-oxide",
        "fileExtension": "svg"
    },
    {
        "ID": 53,
        "file": "Chromium-III.svg",
        "name": "Chromium III",
        "nameSlug": "chromium-iii",
        "fileExtension": "svg"
    },
    {
        "ID": 54,
        "file": "Chromium-II.svg",
        "name": "Chromium II",
        "nameSlug": "chromium-ii",
        "fileExtension": "svg"
    },
    {
        "ID": 55,
        "file": "Chromium-Ion.svg",
        "name": "Chromium Ion",
        "nameSlug": "chromium-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 56,
        "file": "Chromium.svg",
        "name": "Chromium",
        "nameSlug": "chromium",
        "fileExtension": "svg"
    },
    {
        "ID": 57,
        "file": "Copper-III.svg",
        "name": "Copper III",
        "nameSlug": "copper-iii",
        "fileExtension": "svg"
    },
    {
        "ID": 58,
        "file": "Copper-II-Nitrate.svg",
        "name": "Copper II Nitrate",
        "nameSlug": "copper-ii-nitrate",
        "fileExtension": "svg"
    },
    {
        "ID": 59,
        "file": "Copper-II-Sulfate.svg",
        "name": "Copper II Sulfate",
        "nameSlug": "copper-ii-sulfate",
        "fileExtension": "svg"
    },
    {
        "ID": 60,
        "file": "Copper-II.svg",
        "name": "Copper II",
        "nameSlug": "copper-ii",
        "fileExtension": "svg"
    },
    {
        "ID": 61,
        "file": "Copper.svg",
        "name": "Copper",
        "nameSlug": "copper",
        "fileExtension": "svg"
    },
    {
        "ID": 62,
        "file": "Cyanide-Dots.svg",
        "name": "Cyanide Dots",
        "nameSlug": "cyanide-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 63,
        "file": "Cyanide.svg",
        "name": "Cyanide",
        "nameSlug": "cyanide",
        "fileExtension": "svg"
    },
    {
        "ID": 64,
        "file": "DDCF351F.png",
        "name": "DDCF351F",
        "nameSlug": "ddcf351f",
        "fileExtension": "png"
    },
    {
        "ID": 65,
        "file": "Deuterium.svg",
        "name": "Deuterium",
        "nameSlug": "deuterium",
        "fileExtension": "svg"
    },
    {
        "ID": 66,
        "file": "Dinitrogen-Tetroxide.svg",
        "name": "Dinitrogen Tetroxide",
        "nameSlug": "dinitrogen-tetroxide",
        "fileExtension": "svg"
    },
    {
        "ID": 67,
        "file": "Energy.svg",
        "name": "Energy",
        "nameSlug": "energy",
        "fileExtension": "svg"
    },
    {
        "ID": 68,
        "file": "Ethanol.svg",
        "name": "Ethanol",
        "nameSlug": "ethanol",
        "fileExtension": "svg"
    },
    {
        "ID": 69,
        "file": "Ethene.svg",
        "name": "Ethene",
        "nameSlug": "ethene",
        "fileExtension": "svg"
    },
    {
        "ID": 70,
        "file": "filenames.txt",
        "name": "filenames.txt",
        "nameSlug": "filenames.txt",
        "fileExtension": "txt"
    },
    {
        "ID": 71,
        "file": "Fluoride.svg",
        "name": "Fluoride",
        "nameSlug": "fluoride",
        "fileExtension": "svg"
    },
    {
        "ID": 72,
        "file": "Fluorine-Atom.svg",
        "name": "Fluorine Atom",
        "nameSlug": "fluorine-atom",
        "fileExtension": "svg"
    },
    {
        "ID": 73,
        "file": "Fluorine-Ion.svg",
        "name": "Fluorine Ion",
        "nameSlug": "fluorine-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 74,
        "file": "Fluorine.svg",
        "name": "Fluorine",
        "nameSlug": "fluorine",
        "fileExtension": "svg"
    },
    {
        "ID": 75,
        "file": "Generic.svg",
        "name": "Generic",
        "nameSlug": "generic",
        "fileExtension": "svg"
    },
    {
        "ID": 76,
        "file": "Glycerol.svg",
        "name": "Glycerol",
        "nameSlug": "glycerol",
        "fileExtension": "svg"
    },
    {
        "ID": 77,
        "file": "Gold-III.svg",
        "name": "Gold III",
        "nameSlug": "gold-iii",
        "fileExtension": "svg"
    },
    {
        "ID": 78,
        "file": "Gold-I.svg",
        "name": "Gold I",
        "nameSlug": "gold-i",
        "fileExtension": "svg"
    },
    {
        "ID": 79,
        "file": "Gold.svg",
        "name": "Gold",
        "nameSlug": "gold",
        "fileExtension": "svg"
    },
    {
        "ID": 80,
        "file": "Green.svg",
        "name": "Green",
        "nameSlug": "green",
        "fileExtension": "svg"
    },
    {
        "ID": 81,
        "file": "Helium-Ion.svg",
        "name": "Helium Ion",
        "nameSlug": "helium-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 82,
        "file": "Helium.svg",
        "name": "Helium",
        "nameSlug": "helium",
        "fileExtension": "svg"
    },
    {
        "ID": 83,
        "file": "Hydrochloric-Acid-Dots.svg",
        "name": "Hydrochloric Acid Dots",
        "nameSlug": "hydrochloric-acid-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 84,
        "file": "Hydrochloric-Acid.svg",
        "name": "Hydrochloric Acid",
        "nameSlug": "hydrochloric-acid",
        "fileExtension": "svg"
    },
    {
        "ID": 85,
        "file": "Hydrogen-Atom.svg",
        "name": "Hydrogen Atom",
        "nameSlug": "hydrogen-atom",
        "fileExtension": "svg"
    },
    {
        "ID": 86,
        "file": "Hydrogen-Bromide-Dots.svg",
        "name": "Hydrogen Bromide Dots",
        "nameSlug": "hydrogen-bromide-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 87,
        "file": "Hydrogen-Bromide.svg",
        "name": "Hydrogen Bromide",
        "nameSlug": "hydrogen-bromide",
        "fileExtension": "svg"
    },
    {
        "ID": 88,
        "file": "Hydrogen-Cyanide-Dots.svg",
        "name": "Hydrogen Cyanide Dots",
        "nameSlug": "hydrogen-cyanide-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 89,
        "file": "Hydrogen-Cyanide.svg",
        "name": "Hydrogen Cyanide",
        "nameSlug": "hydrogen-cyanide",
        "fileExtension": "svg"
    },
    {
        "ID": 90,
        "file": "Hydrogen-Fluoride.svg",
        "name": "Hydrogen Fluoride",
        "nameSlug": "hydrogen-fluoride",
        "fileExtension": "svg"
    },
    {
        "ID": 91,
        "file": "Hydrogen-Iodide.svg",
        "name": "Hydrogen Iodide",
        "nameSlug": "hydrogen-iodide",
        "fileExtension": "svg"
    },
    {
        "ID": 92,
        "file": "Hydrogen-Ion.svg",
        "name": "Hydrogen Ion",
        "nameSlug": "hydrogen-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 93,
        "file": "Hydrogen-Peroxide.svg",
        "name": "Hydrogen Peroxide",
        "nameSlug": "hydrogen-peroxide",
        "fileExtension": "svg"
    },
    {
        "ID": 94,
        "file": "Hydrogen-Sulfide.svg",
        "name": "Hydrogen Sulfide",
        "nameSlug": "hydrogen-sulfide",
        "fileExtension": "svg"
    },
    {
        "ID": 95,
        "file": "Hydrogen.svg",
        "name": "Hydrogen",
        "nameSlug": "hydrogen",
        "fileExtension": "svg"
    },
    {
        "ID": 96,
        "file": "Hydronium.svg",
        "name": "Hydronium",
        "nameSlug": "hydronium",
        "fileExtension": "svg"
    },
    {
        "ID": 97,
        "file": "Hydroxide-Dots.svg",
        "name": "Hydroxide Dots",
        "nameSlug": "hydroxide-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 98,
        "file": "Hydroxide.svg",
        "name": "Hydroxide",
        "nameSlug": "hydroxide",
        "fileExtension": "svg"
    },
    {
        "ID": 99,
        "file": "Inert.svg",
        "name": "Inert",
        "nameSlug": "inert",
        "fileExtension": "svg"
    },
    {
        "ID": 100,
        "file": "Inhibitor.svg",
        "name": "Inhibitor",
        "nameSlug": "inhibitor",
        "fileExtension": "svg"
    },
    {
        "ID": 101,
        "file": "Iodine-Atom.svg",
        "name": "Iodine Atom",
        "nameSlug": "iodine-atom",
        "fileExtension": "svg"
    },
    {
        "ID": 102,
        "file": "Iodine-Ion.svg",
        "name": "Iodine Ion",
        "nameSlug": "iodine-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 103,
        "file": "Iodine.svg",
        "name": "Iodine",
        "nameSlug": "iodine",
        "fileExtension": "svg"
    },
    {
        "ID": 104,
        "file": "Iron-III.svg",
        "name": "Iron III",
        "nameSlug": "iron-iii",
        "fileExtension": "svg"
    },
    {
        "ID": 105,
        "file": "Iron-II-Sulfate.svg",
        "name": "Iron II Sulfate",
        "nameSlug": "iron-ii-sulfate",
        "fileExtension": "svg"
    },
    {
        "ID": 106,
        "file": "Iron-II.svg",
        "name": "Iron II",
        "nameSlug": "iron-ii",
        "fileExtension": "svg"
    },
    {
        "ID": 107,
        "file": "Iron-I.svg",
        "name": "Iron I",
        "nameSlug": "iron-i",
        "fileExtension": "svg"
    },
    {
        "ID": 108,
        "file": "Iron.svg",
        "name": "Iron",
        "nameSlug": "iron",
        "fileExtension": "svg"
    },
    {
        "ID": 109,
        "file": "Krypton-50.svg",
        "name": "Krypton 50",
        "nameSlug": "krypton-50",
        "fileExtension": "svg"
    },
    {
        "ID": 110,
        "file": "Krypton-89.svg",
        "name": "Krypton 89",
        "nameSlug": "krypton-89",
        "fileExtension": "svg"
    },
    {
        "ID": 111,
        "file": "Lead-Iodide.svg",
        "name": "Lead Iodide",
        "nameSlug": "lead-iodide",
        "fileExtension": "svg"
    },
    {
        "ID": 112,
        "file": "Lead-Ion.svg",
        "name": "Lead Ion",
        "nameSlug": "lead-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 113,
        "file": "Lead.svg",
        "name": "Lead",
        "nameSlug": "lead",
        "fileExtension": "svg"
    },
    {
        "ID": 114,
        "file": "Lithium-Chloride.svg",
        "name": "Lithium Chloride",
        "nameSlug": "lithium-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 115,
        "file": "Lithium-Hydroxide.svg",
        "name": "Lithium Hydroxide",
        "nameSlug": "lithium-hydroxide",
        "fileExtension": "svg"
    },
    {
        "ID": 116,
        "file": "Lithium-Ion.svg",
        "name": "Lithium Ion",
        "nameSlug": "lithium-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 117,
        "file": "Lithium-I.svg",
        "name": "Lithium I",
        "nameSlug": "lithium-i",
        "fileExtension": "svg"
    },
    {
        "ID": 118,
        "file": "Lithium-Nitrate.svg",
        "name": "Lithium Nitrate",
        "nameSlug": "lithium-nitrate",
        "fileExtension": "svg"
    },
    {
        "ID": 119,
        "file": "Lithium-Sulfide.svg",
        "name": "Lithium Sulfide",
        "nameSlug": "lithium-sulfide",
        "fileExtension": "svg"
    },
    {
        "ID": 120,
        "file": "Lithium.svg",
        "name": "Lithium",
        "nameSlug": "lithium",
        "fileExtension": "svg"
    },
    {
        "ID": 121,
        "file": "Magnesium-Ion.svg",
        "name": "Magnesium Ion",
        "nameSlug": "magnesium-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 122,
        "file": "Magnesium-Sulfate.svg",
        "name": "Magnesium Sulfate",
        "nameSlug": "magnesium-sulfate",
        "fileExtension": "svg"
    },
    {
        "ID": 123,
        "file": "Magnesium.svg",
        "name": "Magnesium",
        "nameSlug": "magnesium",
        "fileExtension": "svg"
    },
    {
        "ID": 124,
        "file": "Manganese-Dioxide.svg",
        "name": "Manganese Dioxide",
        "nameSlug": "manganese-dioxide",
        "fileExtension": "svg"
    },
    {
        "ID": 125,
        "file": "Manganese-II-Chloride.svg",
        "name": "Manganese II Chloride",
        "nameSlug": "manganese-ii-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 126,
        "file": "Manganese-II.svg",
        "name": "Manganese II",
        "nameSlug": "manganese-ii",
        "fileExtension": "svg"
    },
    {
        "ID": 127,
        "file": "Manganese-IV.svg",
        "name": "Manganese IV",
        "nameSlug": "manganese-iv",
        "fileExtension": "svg"
    },
    {
        "ID": 128,
        "file": "Manganese.svg",
        "name": "Manganese",
        "nameSlug": "manganese",
        "fileExtension": "svg"
    },
    {
        "ID": 129,
        "file": "Mercury-II.svg",
        "name": "Mercury II",
        "nameSlug": "mercury-ii",
        "fileExtension": "svg"
    },
    {
        "ID": 130,
        "file": "Mercury-I.svg",
        "name": "Mercury I",
        "nameSlug": "mercury-i",
        "fileExtension": "svg"
    },
    {
        "ID": 131,
        "file": "Mercury.svg",
        "name": "Mercury",
        "nameSlug": "mercury",
        "fileExtension": "svg"
    },
    {
        "ID": 132,
        "file": "Methane.svg",
        "name": "Methane",
        "nameSlug": "methane",
        "fileExtension": "svg"
    },
    {
        "ID": 133,
        "file": "Methylamine.svg",
        "name": "Methylamine",
        "nameSlug": "methylamine",
        "fileExtension": "svg"
    },
    {
        "ID": 134,
        "file": "Methylammonium.svg",
        "name": "Methylammonium",
        "nameSlug": "methylammonium",
        "fileExtension": "svg"
    },
    {
        "ID": 135,
        "file": "Neutron-small.svg",
        "name": "Neutron small",
        "nameSlug": "neutron-small",
        "fileExtension": "svg"
    },
    {
        "ID": 136,
        "file": "Neutron.svg",
        "name": "Neutron",
        "nameSlug": "neutron",
        "fileExtension": "svg"
    },
    {
        "ID": 137,
        "file": "Nitrate.svg",
        "name": "Nitrate",
        "nameSlug": "nitrate",
        "fileExtension": "svg"
    },
    {
        "ID": 138,
        "file": "Nitric-Acid.svg",
        "name": "Nitric Acid",
        "nameSlug": "nitric-acid",
        "fileExtension": "svg"
    },
    {
        "ID": 139,
        "file": "Nitric-Oxide.svg",
        "name": "Nitric Oxide",
        "nameSlug": "nitric-oxide",
        "fileExtension": "svg"
    },
    {
        "ID": 140,
        "file": "Nitrogen-Atom.svg",
        "name": "Nitrogen Atom",
        "nameSlug": "nitrogen-atom",
        "fileExtension": "svg"
    },
    {
        "ID": 141,
        "file": "Nitrogen-Dioxide.svg",
        "name": "Nitrogen Dioxide",
        "nameSlug": "nitrogen-dioxide",
        "fileExtension": "svg"
    },
    {
        "ID": 142,
        "file": "Nitrogen-Ion.svg",
        "name": "Nitrogen Ion",
        "nameSlug": "nitrogen-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 143,
        "file": "Nitrogen.svg",
        "name": "Nitrogen",
        "nameSlug": "nitrogen",
        "fileExtension": "svg"
    },
    {
        "ID": 144,
        "file": "Nitrosyl-Chloride.svg",
        "name": "Nitrosyl Chloride",
        "nameSlug": "nitrosyl-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 145,
        "file": "Nitryl-Chloride.svg",
        "name": "Nitryl Chloride",
        "nameSlug": "nitryl-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 146,
        "file": "Nuclear-01.svg",
        "name": "Nuclear 01",
        "nameSlug": "nuclear-01",
        "fileExtension": "svg"
    },
    {
        "ID": 147,
        "file": "Nuclear-02.svg",
        "name": "Nuclear 02",
        "nameSlug": "nuclear-02",
        "fileExtension": "svg"
    },
    {
        "ID": 148,
        "file": "Nuclear-03.svg",
        "name": "Nuclear 03",
        "nameSlug": "nuclear-03",
        "fileExtension": "svg"
    },
    {
        "ID": 149,
        "file": "Nuclear-04.svg",
        "name": "Nuclear 04",
        "nameSlug": "nuclear-04",
        "fileExtension": "svg"
    },
    {
        "ID": 150,
        "file": "Nuclear-05.svg",
        "name": "Nuclear 05",
        "nameSlug": "nuclear-05",
        "fileExtension": "svg"
    },
    {
        "ID": 151,
        "file": "Nuclear-06.svg",
        "name": "Nuclear 06",
        "nameSlug": "nuclear-06",
        "fileExtension": "svg"
    },
    {
        "ID": 152,
        "file": "Nuclear-07.svg",
        "name": "Nuclear 07",
        "nameSlug": "nuclear-07",
        "fileExtension": "svg"
    },
    {
        "ID": 153,
        "file": "Nuclear-p19-01.svg",
        "name": "Nuclear p19 01",
        "nameSlug": "nuclear-p19-01",
        "fileExtension": "svg"
    },
    {
        "ID": 154,
        "file": "Nuclear-p22-01.svg",
        "name": "Nuclear p22 01",
        "nameSlug": "nuclear-p22-01",
        "fileExtension": "svg"
    },
    {
        "ID": 155,
        "file": "Nuclear-p22-05.svg",
        "name": "Nuclear p22 05",
        "nameSlug": "nuclear-p22-05",
        "fileExtension": "svg"
    },
    {
        "ID": 156,
        "file": "Nuclear-p39-01.svg",
        "name": "Nuclear p39 01",
        "nameSlug": "nuclear-p39-01",
        "fileExtension": "svg"
    },
    {
        "ID": 157,
        "file": "Nuclear-p39-02.svg",
        "name": "Nuclear p39 02",
        "nameSlug": "nuclear-p39-02",
        "fileExtension": "svg"
    },
    {
        "ID": 158,
        "file": "Octane.svg",
        "name": "Octane",
        "nameSlug": "octane",
        "fileExtension": "svg"
    },
    {
        "ID": 159,
        "file": "Oxygen-Atom.svg",
        "name": "Oxygen Atom",
        "nameSlug": "oxygen-atom",
        "fileExtension": "svg"
    },
    {
        "ID": 160,
        "file": "Oxygen-Ion.svg",
        "name": "Oxygen Ion",
        "nameSlug": "oxygen-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 161,
        "file": "Oxygen.svg",
        "name": "Oxygen",
        "nameSlug": "oxygen",
        "fileExtension": "svg"
    },
    {
        "ID": 162,
        "file": "Ozone.svg",
        "name": "Ozone",
        "nameSlug": "ozone",
        "fileExtension": "svg"
    },
    {
        "ID": 163,
        "file": "Particle-BetaNegative.svg",
        "name": "Particle BetaNegative",
        "nameSlug": "particle-betanegative",
        "fileExtension": "svg"
    },
    {
        "ID": 164,
        "file": "Pentane.svg",
        "name": "Pentane",
        "nameSlug": "pentane",
        "fileExtension": "svg"
    },
    {
        "ID": 165,
        "file": "Phenylpthalein.svg",
        "name": "Phenylpthalein",
        "nameSlug": "phenylpthalein",
        "fileExtension": "svg"
    },
    {
        "ID": 166,
        "file": "Phosphorus-Ion.svg",
        "name": "Phosphorus Ion",
        "nameSlug": "phosphorus-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 167,
        "file": "Phosphorus-Pentachloride.svg",
        "name": "Phosphorus Pentachloride",
        "nameSlug": "phosphorus-pentachloride",
        "fileExtension": "svg"
    },
    {
        "ID": 168,
        "file": "Phosphorus.svg",
        "name": "Phosphorus",
        "nameSlug": "phosphorus",
        "fileExtension": "svg"
    },
    {
        "ID": 169,
        "file": "Phosphorus-Trichloride.svg",
        "name": "Phosphorus Trichloride",
        "nameSlug": "phosphorus-trichloride",
        "fileExtension": "svg"
    },
    {
        "ID": 170,
        "file": "Potassium-85.svg",
        "name": "Potassium 85",
        "nameSlug": "potassium-85",
        "fileExtension": "svg"
    },
    {
        "ID": 171,
        "file": "Potassium-Bromide.svg",
        "name": "Potassium Bromide",
        "nameSlug": "potassium-bromide",
        "fileExtension": "svg"
    },
    {
        "ID": 172,
        "file": "Potassium-Chloride.svg",
        "name": "Potassium Chloride",
        "nameSlug": "potassium-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 173,
        "file": "Potassium-Ion.svg",
        "name": "Potassium Ion",
        "nameSlug": "potassium-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 174,
        "file": "Potassium-Nitrate.svg",
        "name": "Potassium Nitrate",
        "nameSlug": "potassium-nitrate",
        "fileExtension": "svg"
    },
    {
        "ID": 175,
        "file": "Potassium.svg",
        "name": "Potassium",
        "nameSlug": "potassium",
        "fileExtension": "svg"
    },
    {
        "ID": 176,
        "file": "Propane.svg",
        "name": "Propane",
        "nameSlug": "propane",
        "fileExtension": "svg"
    },
    {
        "ID": 177,
        "file": "Proton-small.svg",
        "name": "Proton small",
        "nameSlug": "proton-small",
        "fileExtension": "svg"
    },
    {
        "ID": 178,
        "file": "Proton.svg",
        "name": "Proton",
        "nameSlug": "proton",
        "fileExtension": "svg"
    },
    {
        "ID": 179,
        "file": "Quark-down-noShadow.svg",
        "name": "Quark down noShadow",
        "nameSlug": "quark-down-noshadow",
        "fileExtension": "svg"
    },
    {
        "ID": 180,
        "file": "Quark-down.svg",
        "name": "Quark down",
        "nameSlug": "quark-down",
        "fileExtension": "svg"
    },
    {
        "ID": 181,
        "file": "Quark-up-noShadow.svg",
        "name": "Quark up noShadow",
        "nameSlug": "quark-up-noshadow",
        "fileExtension": "svg"
    },
    {
        "ID": 182,
        "file": "Quark-up.svg",
        "name": "Quark up",
        "nameSlug": "quark-up",
        "fileExtension": "svg"
    },
    {
        "ID": 183,
        "file": "Red.svg",
        "name": "Red",
        "nameSlug": "red",
        "fileExtension": "svg"
    },
    {
        "ID": 184,
        "file": "Rubidium-96.svg",
        "name": "Rubidium 96",
        "nameSlug": "rubidium-96",
        "fileExtension": "svg"
    },
    {
        "ID": 185,
        "file": "Shape-1.svg",
        "name": "Shape 1",
        "nameSlug": "shape-1",
        "fileExtension": "svg"
    },
    {
        "ID": 186,
        "file": "Shape-2.svg",
        "name": "Shape 2",
        "nameSlug": "shape-2",
        "fileExtension": "svg"
    },
    {
        "ID": 187,
        "file": "Silicon-Dioxide.svg",
        "name": "Silicon Dioxide",
        "nameSlug": "silicon-dioxide",
        "fileExtension": "svg"
    },
    {
        "ID": 188,
        "file": "Silicon-Ion.svg",
        "name": "Silicon Ion",
        "nameSlug": "silicon-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 189,
        "file": "Silicon-IV.svg",
        "name": "Silicon IV",
        "nameSlug": "silicon-iv",
        "fileExtension": "svg"
    },
    {
        "ID": 190,
        "file": "Silicon.svg",
        "name": "Silicon",
        "nameSlug": "silicon",
        "fileExtension": "svg"
    },
    {
        "ID": 191,
        "file": "Silver-Bromide.svg",
        "name": "Silver Bromide",
        "nameSlug": "silver-bromide",
        "fileExtension": "svg"
    },
    {
        "ID": 192,
        "file": "Silver-Carbonate.svg",
        "name": "Silver Carbonate",
        "nameSlug": "silver-carbonate",
        "fileExtension": "svg"
    },
    {
        "ID": 193,
        "file": "Silver-Chloride.svg",
        "name": "Silver Chloride",
        "nameSlug": "silver-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 194,
        "file": "Silver-Hydroxide.svg",
        "name": "Silver Hydroxide",
        "nameSlug": "silver-hydroxide",
        "fileExtension": "svg"
    },
    {
        "ID": 195,
        "file": "Silver-Ion.svg",
        "name": "Silver Ion",
        "nameSlug": "silver-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 196,
        "file": "Silver-Nitrate.svg",
        "name": "Silver Nitrate",
        "nameSlug": "silver-nitrate",
        "fileExtension": "svg"
    },
    {
        "ID": 197,
        "file": "Silver.svg",
        "name": "Silver",
        "nameSlug": "silver",
        "fileExtension": "svg"
    },
    {
        "ID": 198,
        "file": "Sodium-Acetate.svg",
        "name": "Sodium Acetate",
        "nameSlug": "sodium-acetate",
        "fileExtension": "svg"
    },
    {
        "ID": 199,
        "file": "Sodium-Bicarbonate.svg",
        "name": "Sodium Bicarbonate",
        "nameSlug": "sodium-bicarbonate",
        "fileExtension": "svg"
    },
    {
        "ID": 200,
        "file": "Sodium-Carbonate.svg",
        "name": "Sodium Carbonate",
        "nameSlug": "sodium-carbonate",
        "fileExtension": "svg"
    },
    {
        "ID": 201,
        "file": "Sodium-Chloride-Dots.svg",
        "name": "Sodium Chloride Dots",
        "nameSlug": "sodium-chloride-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 202,
        "file": "Sodium-Chloride.svg",
        "name": "Sodium Chloride",
        "nameSlug": "sodium-chloride",
        "fileExtension": "svg"
    },
    {
        "ID": 203,
        "file": "Sodium-Hydroxide-Dots.svg",
        "name": "Sodium Hydroxide Dots",
        "nameSlug": "sodium-hydroxide-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 204,
        "file": "Sodium-Hydroxide.svg",
        "name": "Sodium Hydroxide",
        "nameSlug": "sodium-hydroxide",
        "fileExtension": "svg"
    },
    {
        "ID": 205,
        "file": "Sodium-Ion.svg",
        "name": "Sodium Ion",
        "nameSlug": "sodium-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 206,
        "file": "Sodium-I.svg",
        "name": "Sodium I",
        "nameSlug": "sodium-i",
        "fileExtension": "svg"
    },
    {
        "ID": 207,
        "file": "Sodium-Nitrate.svg",
        "name": "Sodium Nitrate",
        "nameSlug": "sodium-nitrate",
        "fileExtension": "svg"
    },
    {
        "ID": 208,
        "file": "Sodium-Nitrite.svg",
        "name": "Sodium Nitrite",
        "nameSlug": "sodium-nitrite",
        "fileExtension": "svg"
    },
    {
        "ID": 209,
        "file": "Sodium.svg",
        "name": "Sodium",
        "nameSlug": "sodium",
        "fileExtension": "svg"
    },
    {
        "ID": 210,
        "file": "Spark.svg",
        "name": "Spark",
        "nameSlug": "spark",
        "fileExtension": "svg"
    },
    {
        "ID": 211,
        "file": "Sulfate.svg",
        "name": "Sulfate",
        "nameSlug": "sulfate",
        "fileExtension": "svg"
    },
    {
        "ID": 212,
        "file": "Sulfide.svg",
        "name": "Sulfide",
        "nameSlug": "sulfide",
        "fileExtension": "svg"
    },
    {
        "ID": 213,
        "file": "Sulfur-Dioxide.svg",
        "name": "Sulfur Dioxide",
        "nameSlug": "sulfur-dioxide",
        "fileExtension": "svg"
    },
    {
        "ID": 214,
        "file": "Sulfur-Ion.svg",
        "name": "Sulfur Ion",
        "nameSlug": "sulfur-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 215,
        "file": "Sulfur.svg",
        "name": "Sulfur",
        "nameSlug": "sulfur",
        "fileExtension": "svg"
    },
    {
        "ID": 216,
        "file": "Sulfur-Tetrafluoride.svg",
        "name": "Sulfur Tetrafluoride",
        "nameSlug": "sulfur-tetrafluoride",
        "fileExtension": "svg"
    },
    {
        "ID": 217,
        "file": "Thiocyanate.svg",
        "name": "Thiocyanate",
        "nameSlug": "thiocyanate",
        "fileExtension": "svg"
    },
    {
        "ID": 218,
        "file": "Thiocyanic-Acid.svg",
        "name": "Thiocyanic Acid",
        "nameSlug": "thiocyanic-acid",
        "fileExtension": "svg"
    },
    {
        "ID": 219,
        "file": "Tin-II.svg",
        "name": "Tin II",
        "nameSlug": "tin-ii",
        "fileExtension": "svg"
    },
    {
        "ID": 220,
        "file": "Tin-IV.svg",
        "name": "Tin IV",
        "nameSlug": "tin-iv",
        "fileExtension": "svg"
    },
    {
        "ID": 221,
        "file": "Tin.svg",
        "name": "Tin",
        "nameSlug": "tin",
        "fileExtension": "svg"
    },
    {
        "ID": 222,
        "file": "trichlorofluoromethane.svg",
        "name": "Trichlorofluoromethane",
        "nameSlug": "trichlorofluoromethane",
        "fileExtension": "svg"
    },
    {
        "ID": 223,
        "file": "Uranium-235.svg",
        "name": "Uranium-235",
        "nameSlug": "uranium-235",
        "fileExtension": "svg"
    },
    {
        "ID": 224,
        "file": "Water-Dots.svg",
        "name": "Water (Lewis)",
        "nameSlug": "water-dots",
        "fileExtension": "svg"
    },
    {
        "ID": 225,
        "file": "Water.svg",
        "name": "Water",
        "nameSlug": "water",
        "fileExtension": "svg"
    },
    {
        "ID": 226,
        "file": "Zinc-Ion.svg",
        "name": "Zinc Ion",
        "nameSlug": "zinc-ion",
        "fileExtension": "svg"
    },
    {
        "ID": 227,
        "file": "Zinc.svg",
        "name": "Zinc",
        "nameSlug": "zinc",
        "fileExtension": "svg"
    }
]