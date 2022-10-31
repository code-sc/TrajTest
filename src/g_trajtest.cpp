#include <fstream>
#include <gromacs/fileio/oenv.h>
#include <gromacs/fileio/trxio.h>
#include <gromacs/pbcutil/pbc.h>
#include <gromacs/trajectory/trajectoryframe.h>
#include <iostream>
#include <string>

void
positionUpdate (const float dt, rvec* pos, rvec* vel, const rvec* k, int nind)
{
    for (int i = 0; i < nind; i++)
        {
            vel[i][0] = vel[i][0] - k[i][0] * pos[i][0] * dt;
            vel[i][1] = vel[i][1] - k[i][1] * pos[i][1] * dt;
            vel[i][2] = 0;

            pos[i][0] = pos[i][0] + vel[i][0] * dt;
            pos[i][1] = pos[i][1] + vel[i][1] * dt;
            pos[i][2] = pos[i][2] + vel[i][2] * dt;
        }
}

void
writeFrame (std::ofstream& fObj, t_trxframe* fr, t_trxstatus* status, float time, int step,
            rvec* pos, int* ind, int nind)
{
    // Write to dat file
    for (int i = 0; i < nind; i++)
        {
            fObj << time << " " << ind[i] << " " << pos[i][0] << " " << pos[i][1] << " "
                 << pos[i][2] << "\n";
        }

    // Update the frame object and write to XTC
    fr->step = step;
    fr->time = time;
    write_trxframe_indexed (status, fr, nind, ind, nullptr);
    return;
}

void
initFrame (t_trxframe* fr, rvec* pos, int nind)
{
    if (!fr)
        return;

    fr->bStep = TRUE;
    fr->bTime = TRUE;
    fr->bX    = TRUE;
    fr->bPBC  = TRUE;
    fr->bBox  = TRUE;

    fr->step    = 0;
    fr->time    = 0;
    fr->x       = pos;
    fr->pbcType = PbcType::Xyz;

    fr->box[0][0] = 5.0;
    fr->box[1][1] = 5.0;
    fr->box[2][2] = 5.0;

    fr->natoms = nind;
}

void
modelRunner (const float dt, const int steps, int nind)
{
    rvec* pos = nullptr;  // Positions for each atom
    rvec* vel = nullptr;  // Velocity of each atom
    rvec* k   = nullptr;  // Spring constant for each atom
    int* ind  = nullptr;  // Atom indices

    ind = new int[nind];
    pos = new rvec[nind];
    vel = new rvec[nind];
    k   = new rvec[nind];

    // Initialization of positions, velocities and spring constants
    for (int i = 0; i < nind; i++)
        {
            ind[i] = i;

            pos[i][0] = 1.0;
            pos[i][1] = 1.0 + i / 2.0;
            pos[i][2] = 0.0;

            vel[i][0] = 1.0;
            vel[i][1] = 0.0 - 1 / 2.0;
            vel[i][2] = 0.0;

            k[i][0] = 1.0;
            k[i][1] = 1.0;
            k[i][2] = 0.0;
        }

    std::string fName = "test.xtc";

    // XTC setup
    t_trxframe fr;
    clear_trxframe (&fr, TRUE);
    initFrame (&fr, pos, nind);
    t_trxstatus* status = open_trx (fName.c_str (), "w");  // Can only be used to write xtc files

    // Normal file setup
    std::ofstream fObj;
    fObj.open ("trajectory.dat");

    // Simulation loop
    for (int i = 0; i < steps; i++)
        {
            positionUpdate (dt, pos, vel, k, nind);
            writeFrame (fObj, &fr, status, dt * i, i, pos, ind, nind);
        }

    fObj.close ();
    close_trx (status);

    delete[] ind;
    delete[] pos;
    delete[] vel;
    delete[] k;
}

void
readXTC (int nind)
{
    std::ofstream fObj;
    fObj.open ("xtc.dat");
    std::string fName = "test.xtc";

    t_trxframe fr;
    t_trxstatus* status    = nullptr;
    gmx_output_env_t* oenv = nullptr;
    output_env_init_default (&oenv);

    read_first_frame (oenv, &status, fName.c_str (), &fr, TRX_NEED_X);
    /*
       gmx_ouput_env_t* oenv	:
            Pass oenv, where oenv is defined as `gmx_output_env_t* oenv`
            To directly create the oenv object, you'll need `Programcontext` to be passed.
            The output environment is needed to setup the printing of the read frame status
            Data such as he unit of time, the GROMACS header, etc are maintained by oenv.
            With an the default constructor (which explicitly requires `context` ) the the following
            setup is used:

            [1] gmx::IprogramContext& Programcontext = `context` passed to the constructor
            [2] gmx::TimeUnit timeUnit	= gmx::TimeUnit::Picoseconds
            [3] gmx_bool view		= FALSE
            [4] XvgFormat xvgFormat	= XvgFormat::None
            [5] int verbosity		= 0
            [6] trajectory_io_verbosity = 0

            For a proper gromaces program the problem of `context` is solved as the oenv
            initializtion is perforemed by `parse_common_args` which actually also ALOCATES MEMORY
            for the oenv object. `output_env_done` clears the memory allocated by
            `parse_common_args`



       t_trxstatus** status	:
            Pass &status, where status is defined as `t_trxstatus* status = nullptr`
            Will allocate memory during function call

       const char* fn		:
            Pass the file name of the xtc file to open

       t_trxframe* fr		:
            Pass &fr, where fr is defined as `t_trxframe fr`
            Function call will result in `clear_trxframe(fr,TRUE)`

       int flags		:
            use TRX_NEED_X

    */

    for (int i = 0; i < nind; i++)
        {
            fObj << fr.time << " " << i << " " << fr.x[i][0] << " " << fr.x[i][1] << " "
                 << fr.x[i][2] << "\n";
        }

    while (read_next_frame (oenv, status, &fr))
        {
            for (int i = 0; i < nind; i++)
                {
                    fObj << fr.time << " " << i << " " << fr.x[i][0] << " " << fr.x[i][1] << " "
                         << fr.x[i][2] << "\n";
                }
        }

    close_trx (status);  // frees status also
    sfree (fr.x);        // read_first_xtc called by read_first_frame allocates memory for x
    output_env_done (oenv);
    fObj.close ();
    return;
}

int
main (int argc, char* argv[])
{
    float dt  = 0.5;
    int steps = 50000;
    int nind  = 5;

    if (argc == 2)
        {
            if (argv[1][0] == 'w' || argv[1][0] == 'W')
                {
                    modelRunner (dt, steps, nind);
                }
            else if (argv[1][0] == 'r' || argv[1][0] == 'R')
                {
                    readXTC (nind);
                }
            else
                {
                    std::cout << "Invalid options specified. Use\n[1] w/W\twrite "
                                 "trajectory\n[2] r/R\tRead XTC\n";
                    return 1;
                }
        }
    else if (argc == 1)
        {
            modelRunner (dt, steps, nind);
        }
    else
        {
            std::cout << "Invalid number of options specified\n";
            return 1;
        }
    return 0;
}
