# Bachelors Thesis - Design and Optimization of a Low Energy Trajetory to the Moon for Cubesatellites



**Low-Energy Lunar Transfer Trajectory for a 6U CubeSat**

Design and optimization of a low-energy Earth-to-Moon transfer trajectory for a 6U CubeSat equipped with an iodine-fueled micro electric propulsion (EP) system.
The trajectory is formulated as an optimal control problem and solved in two phases using the indirect variational method (Pontryagin's Maximum Principle):

**Stage 1** - Earth to L2: Transfer from low Earth orbit (LEO) to the vicinity of the Earth-Moon L2 libration point. The cost functional minimizes relative selenocentric orbital energy, driving it to approximately −0.12 (non-dimensional) near L2, ensuring passive gravitational capture by the Moon. The optimal thrust arc is localized near L2, with coasting elsewhere.

<img width="505" height="475" alt="Image" src="https://github.com/user-attachments/assets/390ec658-f371-47f0-b643-f8e494d8e864" />

**Stage 2** - Capture & Circularization: Post-capture stabilization (~15 days) uses a heuristic eccentricity-reduction control law to counter rapid orbital evolution driven by lunar and solar perturbations. Subsequent orbit lowering applies a retro-thrust law (thrust anti-parallel to velocity), with engine cutoff enforced above 10,000 km altitude to prevent perigee drop. The spacecraft is progressively wound into a 100 km circular selenocentric orbit.

<img width="485" height="541" alt="Image" src="https://github.com/user-attachments/assets/bd707ccd-31be-4b5c-8f8d-7da148dc0462" />

<img width="489" height="637" alt="Image" src="https://github.com/user-attachments/assets/10aa93b2-b8e8-4e7a-b492-89be85b03ffb" />

<img width="590" height="463" alt="Image" src="https://github.com/user-attachments/assets/55e51e5d-c467-4770-8af5-2bbf9287cc49" />

**Key results:** Total transfer duration ~167 days, propellant mass ~1.8 kg. The method demonstrates feasibility of low-thrust, low-energy lunar insertion for small spacecraft launched on a conventional medium-class vehicle with a kickstage.
