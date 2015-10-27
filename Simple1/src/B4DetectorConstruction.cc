//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B4DetectorConstruction.cc 77601 2013-11-26 17:08:44Z gcosmo $
// 
/// \file B4DetectorConstruction.cc
/// \brief Implementation of the B4DetectorConstruction class

#include "B4DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* B4DetectorConstruction::fMagFieldMessenger = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::B4DetectorConstruction()
 :  G4VUserDetectorConstruction(),
    fAbsorberPV(0),
    fGapPV(0),
    fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::~B4DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::Construct()
{
    // Define materials 
    DefineMaterials();

    // Define volumes
    return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//#######################################
// Define materials
//#######################################

void B4DetectorConstruction::DefineMaterials()
{ 
    // Lead material defined using NIST Manager
    G4NistManager* nistManager = G4NistManager::Instance();
    nistManager->FindOrBuildMaterial("G4_Pb");

    G4double a;  // mass of a mole;
    G4double z;  // z=mean number of protons;  
    G4double density; 

    // From EES    
    G4int ncomponents, natoms;
    G4double fractionmass;
    G4String symbol;


    //#######################################
    // Define elements
    //#######################################

    G4Element* F  = new G4Element("Fluor" ,symbol="F" , z=9.,  a=18.998*g/mole);
    G4Element* Ce = new G4Element("Cerium",symbol="Ce", z=58., a=140.11*g/mole);

    G4Element* B  = new G4Element("Boron"  ,symbol="B" , z=5. , a=10.81*g/mole);
    G4Element* Ge = new G4Element("Germanium" ,symbol="Ge", z=32., a=72.64*g/mole);
    G4Element* O  = new G4Element("Oxygen" ,symbol="O" , z=8. , a=16.00*g/mole);

    G4Element* N  = new G4Element("Nitrogen" ,symbol="N" , z=7., a= 14.01*g/mole);

    G4Element* C  = new G4Element("Carbon" ,symbol="C" , z=6. , a=12.0107*g/mole);
    G4Element* H  = new G4Element("Hydrogen" ,symbol="H" , z=1. , a=1.00794*g/mole);

    G4Element* Si  = new G4Element("Silicon" ,symbol="Si" , z=14. , a=28.085*g/mole);

    G4Element* K = new G4Element("Potassium", symbol="K", z=19. , a=39.098*g/mole);
    G4Element* Cs = new G4Element("Caesium", symbol="Cs", z=55. , a=132.905*g/mole);
    G4Element* Sb  = new G4Element("Antimony" ,symbol="Sb" , z=51. , a=121.760*g/mole);


    //#######################################
    // Define actual materials from elements
    //#######################################

    // Vacuum - from B4 example!
    new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

    // active material
    G4Material* sens = new G4Material("CeF3" , density= 6.16*g/cm3, ncomponents=2);
    sens->AddElement(Ce , natoms=1);
    sens->AddElement(F ,  natoms=3);

    // BGO
    G4Material* bgo = new G4Material("BGO" , density= 7.13*g/cm3, ncomponents=3);
    bgo->AddElement(B   , natoms=4);
    bgo->AddElement(Ge  , natoms=3);
    bgo->AddElement(O   , natoms=12);

    // Tungsten - doesn't need a definition? Copied from EEShashlik
    new G4Material("Tungsten", z=74., a=183.84*g/mole, density=19.25*g/cm3);

    /*
    // From the B4 example
    new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density
    */


    // PMMA fibre material definition
    G4Material* PMMA = new G4Material("PMMA", density=1.18*g/cm3, ncomponents=3);
    PMMA->AddElement(C , natoms=5);
    PMMA->AddElement(O , natoms=2);
    PMMA->AddElement(H , natoms=8);

    // Optical properties - not all meanings clear here, mostly copied
    G4double photonEnergy[] = {
         0.9*eV,  2.00*eV, 2.03*eV, 2.06*eV, 2.09*eV, 2.12*eV, 2.15*eV, 2.18*eV, 2.21*eV, 2.24*eV, 2.27*eV,
         2.30*eV, 2.33*eV, 2.36*eV, 2.39*eV, 2.42*eV, 2.45*eV, 2.48*eV, 2.51*eV, 2.54*eV, 2.57*eV,
         2.60*eV, 2.63*eV, 2.66*eV, 2.69*eV, 2.72*eV, 2.75*eV, 2.78*eV, 2.81*eV, 2.84*eV, 2.87*eV,
         2.90*eV, 2.93*eV, 2.96*eV, 2.99*eV, 3.02*eV, 3.05*eV, 3.08*eV, 3.11*eV, 3.14*eV, 3.17*eV,
         3.20*eV, 3.23*eV, 3.26*eV, 3.29*eV, 3.32*eV, 3.35*eV, 3.38*eV, 3.41*eV, 3.44*eV, 3.47*eV, 5.5*eV
         };

    G4double refractiveIndexWLSfiber[] = {
        1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
        1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
        1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
        1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
        1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49
        };
    assert(sizeof(refractiveIndexWLSfiber) == sizeof(photonEnergy));

    G4double absWLSfiber[] =  {
        20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m,
        20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m,
        20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m,
        20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m,
        20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m, 20.0*m,
        };
    assert(sizeof(absWLSfiber) == sizeof(photonEnergy));




    // Print materials
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//#######################################
// Define volumes
//#######################################

G4VPhysicalVolume* B4DetectorConstruction::DefineVolumes()
{
    // Geometry parameters
    G4int nofLayers = 1;

    /*G4double absoThickness = 10.*mm;
    G4double gapThickness =  5.*mm;
    G4double calorSizeXY  = 10.*cm;*/

    G4double absoThickness = 3.1*mm;
    G4double gapThickness  = 10.*mm;
    G4double calorSizeXY   =2.5*cm;

    // Descriptives
    G4double layerThickness = absoThickness + gapThickness;
    G4double calorThickness = nofLayers * layerThickness;

    // World dimensions
    G4double worldSizeXY = 1.2 * calorSizeXY;
    G4double worldSizeZ  = 1.2 * calorThickness; 

    // Get materials
    G4Material* defaultMaterial = G4Material::GetMaterial("Galactic");
    G4Material* absorberMaterial = G4Material::GetMaterial("Tungsten"); // Was G4_Pb
    G4Material* gapMaterial = G4Material::GetMaterial("CeF3"); // Was liquidArgon

    if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("B4DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
    }  


    //#######################################
    // World
    //#######################################

    G4VSolid* worldS 
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size
                         
    G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
    G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 


    //#######################################
    // Full calorimeter
    // Layers will be placed inside
    //#######################################

    G4VSolid* calorimeterS
    = new G4Box("Calorimeter",     // its name
                 calorSizeXY/2, calorSizeXY/2, calorThickness/2); // its size
                         
    G4LogicalVolume* calorLV
    = new G4LogicalVolume(
                 calorimeterS,     // its solid
                 defaultMaterial,  // its material <-- Essentially an empty block at this point
                 "Calorimeter");   // its name
                                   
    new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 calorLV,          // its logical volume                         
                 "Calorimeter",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 


    //#######################################
    // Layer
    // Constructed from absorber and crystal
    //#######################################

    G4VSolid* layerS 
    = new G4Box("Layer",           // its name
                 calorSizeXY/2, calorSizeXY/2, layerThickness/2); // its size
                         
    G4LogicalVolume* layerLV
    = new G4LogicalVolume(
                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Layer");         // its name

    new G4PVReplica(
                 "Layer",          // its name
                 layerLV,          // its logical volume
                 calorLV,          // its mother
                 kZAxis,           // axis of replication
                 nofLayers,        // number of replica
                 layerThickness);  // witdth of replica


    //=======================================
    // Absorber

    G4VSolid* absorberS 
    = new G4Box("Abso",            // its name
                 calorSizeXY/2, calorSizeXY/2, absoThickness/2); // its size
                         
    G4LogicalVolume* absorberLV
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "Abso");          // its name
                                   
    fAbsorberPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., -gapThickness/2), // its position <-- Does this mean it 'sticks out'?
                 absorberLV,       // its logical volume                         
                 "Abso",           // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 


    //=======================================
    // Gap

    G4VSolid* gapS 
    = new G4Box("Gap",             // its name
                 calorSizeXY/2, calorSizeXY/2, gapThickness/2); // its size
                         
    G4LogicalVolume* gapLV
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "Gap");           // its name
                                   
    fGapPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., absoThickness/2), // its position
                 gapLV,            // its logical volume                         
                 "Gap",            // its name
                 layerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 



    //
    // print parameters
    //
    G4cout << "\n------------------------------------------------------------"
         << "\n---> The calorimeter is " << nofLayers << " layers of: [ "
         << absoThickness/mm << "mm of " << absorberMaterial->GetName() 
         << " + "
         << gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " 
         << "\n------------------------------------------------------------\n";

    //                                        
    // Visualization attributes
    //
    worldLV->SetVisAttributes (G4VisAttributes::Invisible);

    G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(true);

    //calorLV->SetVisAttributes(simpleBoxVisAtt);

    G4VisAttributes* redBox= new G4VisAttributes(G4Colour(1.0,0.0,0.0));
    redBox->SetForceSolid(true);
    absorberLV->SetVisAttributes(redBox);

    G4VisAttributes* yellowBox= new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    yellowBox->SetForceSolid(false);
    gapLV->SetVisAttributes(yellowBox);


    //
    // Always return the physical World
    //
    return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::ConstructSDandField()
{ 
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    G4ThreeVector fieldValue = G4ThreeVector();
    fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
    fMagFieldMessenger->SetVerboseLevel(1);

    // Register the field messenger for deleting
    G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
