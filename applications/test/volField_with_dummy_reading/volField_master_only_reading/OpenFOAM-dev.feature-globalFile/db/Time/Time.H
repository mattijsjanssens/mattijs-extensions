/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::Time

Description
    Class to control time during OpenFOAM simulations that is also the
    top-level objectRegistry.

SourceFiles
    Time.C
    TimeIO.C
    findInstance.C

\*---------------------------------------------------------------------------*/

#ifndef Time_H
#define Time_H

#include "TimePaths.H"
#include "objectRegistry.H"
#include "unwatchedIOdictionary.H"
#include "FIFOStack.H"
#include "clock.H"
#include "cpuTime.H"
#include "TimeState.H"
#include "Switch.H"
#include "instantList.H"
#include "NamedEnum.H"
#include "typeInfo.H"
#include "dlLibraryTable.H"
#include "functionObjectList.H"
#include "fileMonitor.H"
#include "sigWriteNow.H"
#include "sigStopAtWriteNow.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class argList;

/*---------------------------------------------------------------------------*\
                             Class Time Declaration
\*---------------------------------------------------------------------------*/

class Time
:
    public clock,
    public cpuTime,
    public TimePaths,
    public objectRegistry,
    public TimeState
{
    // Private data

        //- file-change monitor for all registered files
        mutable autoPtr<fileMonitor> monitorPtr_;

        //- Any loaded dynamic libraries. Make sure to construct before
        //  reading controlDict.
        dlLibraryTable libs_;

        //- The controlDict
        unwatchedIOdictionary controlDict_;

public:

        //- Write control options
        enum writeControls
        {
            wcTimeStep,
            wcRunTime,
            wcAdjustableRunTime,
            wcClockTime,
            wcCpuTime
        };

        //- Stop-run control options
        enum stopAtControls
        {
            saEndTime,    //!< stop when Time reaches the prescribed endTime
            saNoWriteNow, //!< set endTime to stop immediately w/o writing
            saWriteNow,   //!< set endTime to stop immediately w/ writing
            saNextWrite   //!< stop the next time data are written
        };

        //- Supported time directory name formats
        enum fmtflags
        {
            general    = 0,
            fixed      = ios_base::fixed,
            scientific = ios_base::scientific
        };


protected:

    // Protected data

        label  startTimeIndex_;
        scalar startTime_;
        mutable scalar endTime_;

        static const NamedEnum<stopAtControls, 4> stopAtControlNames_;
        mutable stopAtControls stopAt_;

        static const NamedEnum<writeControls, 5> writeControlNames_;
        writeControls writeControl_;

        scalar writeInterval_;

        label  purgeWrite_;
        mutable FIFOStack<word> previousWriteTimes_;

        // One-shot writing
        bool writeOnce_;

        //- Is the time currently being sub-cycled?
        bool subCycling_;

        //- If time is being sub-cycled this is the previous TimeState
        autoPtr<TimeState> prevTimeState_;


        // Signal handlers for writing

            //- Enable one-shot writing upon signal
            sigWriteNow sigWriteNow_;

            //- Enable write and clean exit upon signal
            sigStopAtWriteNow sigStopAtWriteNow_;


        //- Time directory name format
        static fmtflags format_;

        //- Time directory name precision
        static int precision_;

        //- Maximum time directory name precision
        static const int maxPrecision_;

        //- Adjust the time step so that writing occurs at the specified time
        void adjustDeltaT();

        //- Set the controls from the current controlDict
        void setControls();

        //- Read the control dictionary and set the write controls etc.
        virtual void readDict();


private:

        //- Default write option
        IOstream::streamFormat writeFormat_;

        //- Default output file format version number
        IOstream::versionNumber writeVersion_;

        //- Default output compression
        IOstream::compressionType writeCompression_;

        //- Default graph format
        word graphFormat_;

        //- Is runtime modification of dictionaries allowed?
        Switch runTimeModifiable_;

        //- Function objects executed at start and on ++, +=
        mutable functionObjectList functionObjects_;


public:

    TypeName("time");

    //- The default control dictionary name (normally "controlDict")
    static word controlDictName;


    // Constructors

        //- Construct given name of dictionary to read and argument list
        Time
        (
            const word& name,
            const argList& args,
            const word& systemName = "system",
            const word& constantName = "constant"
        );

        //- Construct given name of dictionary to read, rootPath and casePath
        Time
        (
            const word& name,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant",
            const bool enableFunctionObjects = true
        );

        //- Construct given dictionary, rootPath and casePath
        Time
        (
            const dictionary& dict,
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant",
            const bool enableFunctionObjects = true
        );

        //- Construct given endTime, rootPath and casePath
        Time
        (
            const fileName& rootPath,
            const fileName& caseName,
            const word& systemName = "system",
            const word& constantName = "constant",
            const bool enableFunctionObjects = true
        );


    //- Destructor
    virtual ~Time();


    // Member functions

        // Database functions

            //- Return root path
            const fileName& rootPath() const
            {
                return TimePaths::rootPath();
            }

            //- Return case name
            const fileName& caseName() const
            {
                return TimePaths::caseName();
            }

            //- Return path
            fileName path() const
            {
                return rootPath()/caseName();
            }

            const dictionary& controlDict() const
            {
                return controlDict_;
            }

            virtual const fileName& dbDir() const
            {
                return fileName::null;
            }

            //- Return current time path
            fileName timePath() const
            {
                return path()/timeName();
            }

            //- Default write format
            IOstream::streamFormat writeFormat() const
            {
                return writeFormat_;
            }

            //- Default write version number
            IOstream::versionNumber writeVersion() const
            {
                return writeVersion_;
            }

            //- Default write compression
            IOstream::compressionType writeCompression() const
            {
                return writeCompression_;
            }

            //- Default graph format
            const word& graphFormat() const
            {
                return graphFormat_;
            }

            //- Supports re-reading
            const Switch& runTimeModifiable() const
            {
                return runTimeModifiable_;
            }

            //- Read control dictionary, update controls and time
            virtual bool read();

            // Automatic rereading

                //- Read the objects that have been modified
                void readModifiedObjects();

                //- Helper: add watches for list of files
                void addWatches(regIOobject&, const fileNameList&) const;

                //- Find index (or -1) of file in list of handles
                label findWatch
                (
                    const labelList& watchIndices,
                    const fileName&
                ) const;

                //- Add watching of a file. Returns handle
                label addTimeWatch(const fileName&) const;

                //- Remove watch on a file (using handle)
                bool removeWatch(const label) const;

                //- Get name of file being watched (using handle)
                const fileName& getFile(const label) const;

                //- Get current state of file (using handle)
                fileMonitor::fileState getState(const label) const;

                //- Set current state of file (using handle) to unmodified
                void setUnmodified(const label) const;


            //- Return the location of "dir" containing the file "name".
            //  (eg, used in reading mesh data)
            //  If name is null, search for the directory "dir" only.
            //  Does not search beyond stopInstance (if set) or constant.
            word findInstance
            (
                const fileName& dir,
                const word& name = word::null,
                const IOobject::readOption rOpt = IOobject::MUST_READ,
                const word& stopInstance = word::null
            ) const;

            //- Search the case for valid time directories
            instantList times() const;

            //- Search the case for the time directory path
            //  corresponding to the given instance
            word findInstancePath(const fileName& path, const instant&) const;

            //- Search the case for the time directory path
            //  corresponding to the given instance
            word findInstancePath(const instant&) const;

            //- Search the case for the time closest to the given time
            instant findClosestTime(const scalar) const;

            //- Search instantList for the time index closest to the given time
            static label findClosestTimeIndex
            (
                const instantList&,
                const scalar,
                const word& constantName = "constant"
            );

            //- Write time dictionary to the \<time\>/uniform directory
            virtual bool writeTimeDict() const;

            //- Write using given format, version and compression
            virtual bool writeObject
            (
                IOstream::streamFormat,
                IOstream::versionNumber,
                IOstream::compressionType
            ) const;

            //- Write the objects now (not at end of iteration) and continue
            //  the run
            bool writeNow();

            //- Write the objects now (not at end of iteration) and end the run
            bool writeAndEnd();

            //- Write the objects once (one shot) and continue the run
            void writeOnce();


        // Access

            //- Return time name of given scalar time
            //  formatted with given precision
            static word timeName
            (
                const scalar,
                const int precision = precision_
            );

            //- Return current time name
            virtual word timeName() const;

            //- Search a given directory for valid time directories
            static instantList findTimes
            (
                const fileName&,
                const word& constantName = "constant"
            );

            //- Return start time index
            virtual label startTimeIndex() const;

            //- Return start time
            virtual dimensionedScalar startTime() const;

            //- Return end time
            virtual dimensionedScalar endTime() const;

            //- Return the list of function objects
            const functionObjectList& functionObjects() const
            {
                return functionObjects_;
            }

            //- External access to the loaded libraries
            const dlLibraryTable& libs() const
            {
                return libs_;
            }

            //- External access to the loaded libraries
            dlLibraryTable& libs()
            {
                return libs_;
            }

            //- Return true if time currently being sub-cycled, otherwise false
            bool subCycling() const
            {
                return subCycling_;
            }

            //- Return previous TimeState if time is being sub-cycled
            const TimeState& prevTimeState() const
            {
                return prevTimeState_();
            }


        // Check

            //- Return true if run should continue,
            //  also invokes the functionObjectList::end() method
            //  when the time goes out of range
            //  \note
            //  For correct behaviour, the following style of time-loop
            //  is recommended:
            //  \code
            //      while (runTime.run())
            //      {
            //          runTime++;
            //          solve;
            //          runTime.write();
            //      }
            //  \endcode
            virtual bool run() const;

            //- Return true if run should continue and if so increment time
            //  also invokes the functionObjectList::end() method
            //  when the time goes out of range
            //  \note
            //  For correct behaviour, the following style of time-loop
            //  is recommended:
            //  \code
            //      while (runTime.loop())
            //      {
            //          solve;
            //          runTime.write();
            //      }
            //  \endcode
            virtual bool loop();

            //- Return true if end of run,
            //  does not invoke any functionObject methods
            //  \note
            //      The rounding heuristics near endTime mean that
            //      \code run() \endcode and \code !end() \endcode may
            //      not yield the same result
            virtual bool end() const;


        // Edit

            //- Adjust the current stopAtControl. Note that this value
            //  only persists until the next time the dictionary is read.
            //  Return true if the stopAtControl changed.
            virtual bool stopAt(const stopAtControls) const;

            //- Reset the time and time-index to those of the given time
            virtual void setTime(const Time&);

            //- Reset the time and time-index
            virtual void setTime(const instant&, const label newIndex);

            //- Reset the time and time-index
            virtual void setTime
            (
                const dimensionedScalar&,
                const label newIndex
            );

            //- Reset the time and time-index
            virtual void setTime(const scalar, const label newIndex);

            //- Reset end time
            virtual void setEndTime(const dimensionedScalar&);

            //- Reset end time
            virtual void setEndTime(const scalar);

            //- Reset time step
            virtual void setDeltaT
            (
                const dimensionedScalar&,
                const bool adjustDeltaT = true
            );

            //- Reset time step
            virtual void setDeltaT
            (
                const scalar,
                const bool adjustDeltaT = true
            );

            //- Set time to sub-cycle for the given number of steps
            virtual TimeState subCycle(const label nSubCycles);

            //- Reset time after sub-cycling back to previous TimeState
            virtual void endSubCycle();

            //- Return non-const access to the list of function objects
            functionObjectList& functionObjects()
            {
                return functionObjects_;
            }


    // Member operators

        //- Set deltaT to that specified and increment time via operator++()
        virtual Time& operator+=(const dimensionedScalar&);

        //- Set deltaT to that specified and increment time via operator++()
        virtual Time& operator+=(const scalar);

        //- Prefix increment,
        //  also invokes the functionObjectList::start() or
        //  functionObjectList::execute() method, depending on the time-index
        virtual Time& operator++();

        //- Postfix increment, this is identical to the prefix increment
        virtual Time& operator++(int);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
