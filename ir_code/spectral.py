import random


# ============================================================
# HIGH-LEVEL PLAYER DATABASE
# ============================================================

HIGH_LEVEL_PLAYERS = {
    "Stephen Afolabi 2",
    "Emmanuel Okoli",
    "Lembo Muize",
    "Brian",
    "papi",
    "Chris Kalubo",
    "Daryl Adopo",
    "The Joker",
    "Bobby Ralph",
    "Fawas",
    "Gr8tDane",
    "KC O",
    "Leef",
}


POSITIONS = [
    "Defender",
    "Midfielder",
    "Forward"
]


# ============================================================
# INPUT VALIDATION
# ============================================================

def get_positive_integer(prompt, allow_zero=False):

    while True:

        try:

            value = int(input(prompt).strip())

            if allow_zero:

                if value < 0:

                    print()
                    print(
                        "❌ Please enter 0 or a positive number."
                    )

                    continue

            else:

                if value <= 0:

                    print()
                    print(
                        "❌ Please enter a number greater than 0."
                    )

                    continue

            return value

        except ValueError:

            print()
            print("❌ Invalid input.")
            print("Please enter a whole number.")

        except KeyboardInterrupt:

            print("\n\nProgram cancelled by user.")
            raise SystemExit


# ============================================================
# HIGH-LEVEL PLAYER CHECK
# ============================================================

def is_high_level_player(name):

    """
    Determine whether a player is in the high-level
    player database.

    Comparison is:
        - case-insensitive
        - whitespace-insensitive at the beginning/end
    """

    name_clean = name.strip().casefold()

    return any(
        player.strip().casefold() == name_clean
        for player in HIGH_LEVEL_PLAYERS
    )


# ============================================================
# PLAYER INPUT
# ============================================================

def get_player_name(position, number, existing_players):

    while True:

        try:

            name = input(
                f"Enter the name of {position.lower()} {number}: "
            ).strip()

        except KeyboardInterrupt:

            print("\n\nProgram cancelled by user.")
            raise SystemExit

        # ----------------------------------------------------
        # EMPTY NAME
        # ----------------------------------------------------

        if not name:

            print()
            print(
                "❌ Player name cannot be empty."
            )

            continue

        # ----------------------------------------------------
        # NAME LENGTH
        # ----------------------------------------------------

        if len(name) > 50:

            print()
            print(
                "❌ Player name is too long."
            )

            print(
                "Please enter a name with 50 characters "
                "or fewer."
            )

            continue

        # ----------------------------------------------------
        # DUPLICATE CHECK
        # ----------------------------------------------------

        duplicate = any(
            player["name"].casefold() == name.casefold()
            for player in existing_players
        )

        if duplicate:

            print()
            print(
                f"❌ '{name}' has already been entered."
            )

            print(
                "Please enter a different player name."
            )

            continue

        # ----------------------------------------------------
        # HIGH-LEVEL CLASSIFICATION
        # ----------------------------------------------------

        if is_high_level_player(name):

            print()
            print(
                f"⭐ {name} has been identified as a "
                f"HIGH-LEVEL PLAYER."
            )

        else:

            print()
            print(
                f"✓ {name} has been identified as a "
                f"regular player."
            )

        return name


# ============================================================
# CREATE PLAYER
# ============================================================

def create_player(name, position):

    high_level = is_high_level_player(name)

    return {
        "name": name,
        "position": position,
        "high_level": high_level
    }


# ============================================================
# SEPARATE PLAYERS
# ============================================================

def separate_players(players):

    high_level_players = []
    regular_players = []

    for player in players:

        if player["high_level"]:

            high_level_players.append(player)

        else:

            regular_players.append(player)

    return high_level_players, regular_players


# ============================================================
# GROUP PLAYERS BY POSITION
# ============================================================

def group_by_position(players):

    grouped = {
        "Defender": [],
        "Midfielder": [],
        "Forward": []
    }

    for player in players:

        position = player["position"]

        if position in grouped:

            grouped[position].append(player)

    return grouped


# ============================================================
# COUNT PLAYERS BY POSITION
# ============================================================

def count_position(team, position):

    return sum(
        1
        for player in team
        if player["position"] == position
    )


# ============================================================
# COUNT HIGH-LEVEL PLAYERS
# ============================================================

def count_high_level_players(team):

    return sum(
        1
        for player in team
        if player["high_level"]
    )


# ============================================================
# COUNT HIGH-LEVEL PLAYERS BY POSITION
# ============================================================

def count_high_level_position(team, position):

    return sum(
        1
        for player in team
        if (
            player["position"] == position
            and player["high_level"]
        )
    )


# ============================================================
# DISTRIBUTE HIGH-LEVEL PLAYERS
# ============================================================

def distribute_high_level_players(
    high_level_players,
    black_team,
    white_team
):

    """
    Distribute HIGH-LEVEL PLAYERS independently by position.

    This is the most important part of the algorithm.

    Example:

        4 high-level defenders
            Black = 2
            White = 2

        3 high-level midfielders
            Black = 2
            White = 1

        2 high-level forwards
            Black = 1
            White = 1

    Therefore, high-level players of one position cannot
    become concentrated on a single team.

    The difference between Black and White for any position
    will never be greater than 1.
    """

    players_by_position = group_by_position(
        high_level_players
    )

    # --------------------------------------------------------
    # RANDOMIZE EACH POSITION INDEPENDENTLY
    # --------------------------------------------------------

    for position in POSITIONS:

        random.shuffle(
            players_by_position[position]
        )

    # --------------------------------------------------------
    # DISTRIBUTE POSITION BY POSITION
    # --------------------------------------------------------

    for position in POSITIONS:

        players = players_by_position[position]

        for player in players:

            black_count = count_high_level_position(
                black_team,
                position
            )

            white_count = count_high_level_position(
                white_team,
                position
            )

            # ------------------------------------------------
            # Give player to team with fewer high-level
            # players at this particular position.
            # ------------------------------------------------

            if black_count < white_count:

                black_team.append(player)

            elif white_count < black_count:

                white_team.append(player)

            else:

                # ------------------------------------------------
                # Equal → random decision.
                # ------------------------------------------------

                if random.choice([True, False]):

                    black_team.append(player)

                else:

                    white_team.append(player)


# ============================================================
# CALCULATE TARGET POSITION DISTRIBUTION
# ============================================================

def calculate_target_position_counts(
    all_players,
    total_black_target
):

    """
    Determine how many players of each position should
    ideally be assigned to Black.

    Each position is kept as balanced as possible while
    ensuring the overall Black team reaches its target size.

    Example:

        Total defenders   = 5
        Total midfielders = 6
        Total forwards    = 5

        Total players = 16
        Black target = 8

    Possible target:

        Defenders   = 3
        Midfielders = 3
        Forwards    = 2

        Total = 8
    """

    players_by_position = group_by_position(
        all_players
    )

    targets = {}

    # --------------------------------------------------------
    # Start with floor(count / 2)
    # --------------------------------------------------------

    for position in POSITIONS:

        total_position_players = len(
            players_by_position[position]
        )

        targets[position] = (
            total_position_players // 2
        )

    current_total = sum(
        targets.values()
    )

    additional_players_needed = (
        total_black_target - current_total
    )

    # --------------------------------------------------------
    # Positions that have an odd number of players can
    # receive the extra player.
    # --------------------------------------------------------

    odd_positions = [
        position
        for position in POSITIONS
        if len(players_by_position[position]) % 2 == 1
    ]

    random.shuffle(odd_positions)

    # --------------------------------------------------------
    # Add extra players until Black reaches target size.
    # --------------------------------------------------------

    for position in odd_positions:

        if additional_players_needed <= 0:

            break

        targets[position] += 1

        additional_players_needed -= 1

    return targets


# ============================================================
# DISTRIBUTE REGULAR PLAYERS
# ============================================================

def distribute_regular_players(
    regular_players,
    black_team,
    white_team,
    all_players
):

    """
    Distribute regular players while maintaining:

        1. Overall team size
        2. Position balance
        3. Previously established high-level distribution

    High-level players are distributed first.
    Regular players then fill the remaining positions.
    """

    total_players = len(all_players)

    black_target_size = (
        (total_players + 1) // 2
    )

    white_target_size = (
        total_players // 2
    )

    # --------------------------------------------------------
    # Determine target number of each position for Black.
    # --------------------------------------------------------

    black_position_targets = (
        calculate_target_position_counts(
            all_players,
            black_target_size
        )
    )

    # --------------------------------------------------------
    # Group regular players by position.
    # --------------------------------------------------------

    regular_by_position = group_by_position(
        regular_players
    )

    # --------------------------------------------------------
    # Randomize regular players.
    # --------------------------------------------------------

    for position in POSITIONS:

        random.shuffle(
            regular_by_position[position]
        )

    # --------------------------------------------------------
    # Fill each position.
    # --------------------------------------------------------

    for position in POSITIONS:

        players = regular_by_position[position]

        for player in players:

            black_position_count = count_position(
                black_team,
                position
            )

            white_position_count = count_position(
                white_team,
                position
            )

            black_total = len(black_team)
            white_total = len(white_team)

            black_position_target = (
                black_position_targets[position]
            )

            # ------------------------------------------------
            # If Black still needs this position, assign there.
            # ------------------------------------------------

            if (
                black_position_count
                < black_position_target
            ):

                black_team.append(player)

                continue

            # ------------------------------------------------
            # Otherwise assign to White.
            # ------------------------------------------------

            if white_total < white_target_size:

                white_team.append(player)

                continue

            # ------------------------------------------------
            # Safety fallback
            # ------------------------------------------------

            if black_total < black_target_size:

                black_team.append(player)

            else:

                white_team.append(player)


# ============================================================
# DISPLAY TEAM
# ============================================================

def display_team(team, team_name):

    print()
    print("=" * 80)
    print(team_name)
    print("=" * 80)

    if not team:

        print("No players assigned.")

        return

    for number, player in enumerate(
        team,
        start=1
    ):

        level = (
            "⭐ HIGH-LEVEL"
            if player["high_level"]
            else "Regular"
        )

        print(
            f"{number:2}. "
            f"{player['name']:<25} "
            f"{player['position']:<15} "
            f"{level}"
        )


# ============================================================
# DISPLAY HIGH-LEVEL SUMMARY
# ============================================================

def display_high_level_summary(
    all_players,
    black_team,
    white_team
):

    high_level_players = [
        player
        for player in all_players
        if player["high_level"]
    ]

    print()
    print("=" * 80)
    print("HIGH-LEVEL PLAYER IDENTIFICATION")
    print("=" * 80)

    if not high_level_players:

        print(
            "No entered players were found in the "
            "high-level player database."
        )

        return

    print(
        f"High-level players identified: "
        f"{len(high_level_players)}"
    )

    print()

    for number, player in enumerate(
        high_level_players,
        start=1
    ):

        if player in black_team:

            assigned_team = "BLACK TEAM"

        elif player in white_team:

            assigned_team = "WHITE TEAM"

        else:

            assigned_team = "UNASSIGNED"

        print(
            f"{number:2}. "
            f"{player['name']:<25} "
            f"{player['position']:<15} "
            f"→ {assigned_team}"
        )


# ============================================================
# DISPLAY HIGH-LEVEL POSITION BALANCE
# ============================================================

def display_high_level_position_balance(
    black_team,
    white_team
):

    print()
    print("=" * 80)
    print("HIGH-LEVEL PLAYER POSITION BALANCE")
    print("=" * 80)

    print()
    print(
        f"{'POSITION':<15}"
        f"{'BLACK':<15}"
        f"{'WHITE':<15}"
        f"{'DIFFERENCE':<15}"
    )

    print("-" * 60)

    for position in POSITIONS:

        black_count = count_high_level_position(
            black_team,
            position
        )

        white_count = count_high_level_position(
            white_team,
            position
        )

        difference = abs(
            black_count - white_count
        )

        print(
            f"{position:<15}"
            f"{black_count:<15}"
            f"{white_count:<15}"
            f"{difference:<15}"
        )

    print()

    perfect = True

    for position in POSITIONS:

        black_count = count_high_level_position(
            black_team,
            position
        )

        white_count = count_high_level_position(
            white_team,
            position
        )

        if abs(
            black_count - white_count
        ) > 1:

            perfect = False

    if perfect:

        print(
            "✅ High-level players are balanced by "
            "position."
        )

        print(
            "✅ No position has more than a one-player "
            "high-level difference."
        )

    else:

        print(
            "⚠️ High-level players could not be "
            "perfectly balanced."
        )


# ============================================================
# DISPLAY POSITION SUMMARY
# ============================================================

def display_position_summary(
    black_team,
    white_team
):

    print()
    print("=" * 80)
    print("OVERALL POSITION SUMMARY")
    print("=" * 80)

    print()
    print(
        f"{'POSITION':<15}"
        f"{'BLACK':<15}"
        f"{'WHITE':<15}"
        f"{'DIFFERENCE':<15}"
    )

    print("-" * 60)

    for position in POSITIONS:

        black_count = count_position(
            black_team,
            position
        )

        white_count = count_position(
            white_team,
            position
        )

        difference = abs(
            black_count - white_count
        )

        print(
            f"{position:<15}"
            f"{black_count:<15}"
            f"{white_count:<15}"
            f"{difference:<15}"
        )


# ============================================================
# DISPLAY BALANCE CHECK
# ============================================================

def display_balance_check(
    black_team,
    white_team
):

    team_difference = abs(
        len(black_team)
        - len(white_team)
    )

    black_high_level = count_high_level_players(
        black_team
    )

    white_high_level = count_high_level_players(
        white_team
    )

    high_level_difference = abs(
        black_high_level
        - white_high_level
    )

    print()
    print("=" * 80)
    print("TEAM BALANCE CHECK")
    print("=" * 80)

    # --------------------------------------------------------
    # TEAM SIZE
    # --------------------------------------------------------

    print()
    print(
        f"Black Team : {len(black_team)} players"
    )

    print(
        f"White Team : {len(white_team)} players"
    )

    print(
        f"Player difference: {team_difference}"
    )

    # --------------------------------------------------------
    # HIGH-LEVEL TOTAL
    # --------------------------------------------------------

    print()

    print(
        f"Black high-level players : "
        f"{black_high_level}"
    )

    print(
        f"White high-level players : "
        f"{white_high_level}"
    )

    print(
        f"High-level difference: "
        f"{high_level_difference}"
    )

    # --------------------------------------------------------
    # POSITION BALANCE
    # --------------------------------------------------------

    print()
    print("POSITION BALANCE")

    position_balanced = True

    for position in POSITIONS:

        black_count = count_position(
            black_team,
            position
        )

        white_count = count_position(
            white_team,
            position
        )

        difference = abs(
            black_count - white_count
        )

        print(
            f"{position:<12}: "
            f"Black {black_count} | "
            f"White {white_count} | "
            f"Difference {difference}"
        )

        if difference > 1:

            position_balanced = False

    # --------------------------------------------------------
    # FINAL ASSESSMENT
    # --------------------------------------------------------

    print()
    print("-" * 80)

    if team_difference == 0:

        print(
            "✅ Teams have exactly the same number "
            "of players."
        )

    elif team_difference == 1:

        print(
            "✅ Teams are balanced as closely as "
            "possible by total player count."
        )

    else:

        print(
            "⚠️ Teams have more than one player "
            "difference."
        )

    if high_level_difference == 0:

        print(
            "⭐ HIGH-LEVEL PLAYERS ARE PERFECTLY "
            "BALANCED OVERALL."
        )

    elif high_level_difference == 1:

        print(
            "⭐ HIGH-LEVEL PLAYERS ARE BALANCED "
            "AS CLOSELY AS POSSIBLE OVERALL."
        )

    else:

        print(
            "⚠️ High-level players could not be "
            "balanced equally overall."
        )

    if position_balanced:

        print(
            "⚽ OVERALL POSITIONS ARE BALANCED "
            "AS CLOSELY AS POSSIBLE."
        )

    else:

        print(
            "⚠️ Some positions have more than one "
            "player difference."
        )


# ============================================================
# DISPLAY COMPLETE HIGH-LEVEL POSITION TABLE
# ============================================================

def display_high_level_detail(
    black_team,
    white_team
):

    print()
    print("=" * 80)
    print("HIGH-LEVEL PLAYERS BY POSITION")
    print("=" * 80)

    for position in POSITIONS:

        black_players = [
            player["name"]
            for player in black_team
            if (
                player["position"] == position
                and player["high_level"]
            )
        ]

        white_players = [
            player["name"]
            for player in white_team
            if (
                player["position"] == position
                and player["high_level"]
            )
        ]

        print()
        print(f"{position.upper()}")

        print(
            f"  BLACK ({len(black_players)}):"
        )

        if black_players:

            for name in black_players:

                print(
                    f"      ⭐ {name}"
                )

        else:

            print(
                "      None"
            )

        print(
            f"  WHITE ({len(white_players)}):"
        )

        if white_players:

            for name in white_players:

                print(
                    f"      ⭐ {name}"
                )

        else:

            print(
                "      None"
            )


# ============================================================
# MAIN PLAYER GROUPING FUNCTION
# ============================================================

def Players_Grouping():

    try:

        defenders = []
        midfielders = []
        forwards = []

        all_players = []

        Black_team = []
        White_team = []

        # ====================================================
        # PROGRAM HEADER
        # ====================================================

        print()
        print("=" * 80)
        print(
            "             RANDOM FOOTBALL TEAM GROUPING"
        )
        print("=" * 80)

        print()

        print(
            "The system will automatically identify "
            "high-level players."
        )

        print(
            f"High-level player database contains "
            f"{len(HIGH_LEVEL_PLAYERS)} player(s)."
        )

        # ====================================================
        # TOTAL PLAYERS
        # ====================================================

        players_num = get_positive_integer(
            "Enter the total number of players: "
        )

        # ====================================================
        # DEFENDERS
        # ====================================================

        print()
        print("-" * 80)
        print("DEFENDERS")
        print("-" * 80)

        while True:

            players_defenders = get_positive_integer(
                "Enter the number of defenders: ",
                allow_zero=True
            )

            if players_defenders > players_num:

                print()

                print(
                    "❌ You cannot have more defenders "
                    "than total players."
                )

                continue

            break

        for i in range(players_defenders):

            name = get_player_name(
                "defender",
                i + 1,
                all_players
            )

            player = create_player(
                name,
                "Defender"
            )

            defenders.append(player)
            all_players.append(player)

        # ====================================================
        # MIDFIELDERS
        # ====================================================

        remaining = (
            players_num
            - players_defenders
        )

        print()
        print("-" * 80)
        print("MIDFIELDERS")
        print("-" * 80)

        while True:

            players_midfielders = get_positive_integer(
                f"Enter the number of midfielders "
                f"(0-{remaining}): ",
                allow_zero=True
            )

            if players_midfielders > remaining:

                print()

                print(
                    f"❌ You only have {remaining} "
                    f"player(s) remaining."
                )

                continue

            break

        for i in range(players_midfielders):

            name = get_player_name(
                "midfielder",
                i + 1,
                all_players
            )

            player = create_player(
                name,
                "Midfielder"
            )

            midfielders.append(player)
            all_players.append(player)

        # ====================================================
        # FORWARDS
        # ====================================================

        remaining = (
            players_num
            - players_defenders
            - players_midfielders
        )

        print()
        print("-" * 80)
        print("FORWARDS")
        print("-" * 80)

        while True:

            players_forwards = get_positive_integer(
                f"Enter the number of forwards "
                f"(0-{remaining}): ",
                allow_zero=True
            )

            if players_forwards > remaining:

                print()

                print(
                    f"❌ You only have {remaining} "
                    f"player(s) remaining."
                )

                continue

            break

        for i in range(players_forwards):

            name = get_player_name(
                "forward",
                i + 1,
                all_players
            )

            player = create_player(
                name,
                "Forward"
            )

            forwards.append(player)
            all_players.append(player)

        # ====================================================
        # VERIFY PLAYER COUNT
        # ====================================================

        entered_players = len(all_players)

        print()
        print("=" * 80)
        print("PLAYER COUNT")
        print("=" * 80)

        print(
            f"Expected players : {players_num}"
        )

        print(
            f"Entered players  : {entered_players}"
        )

        if entered_players != players_num:

            print()

            print(
                "❌ ERROR: Player count does not match."
            )

            print(
                "The teams cannot be created."
            )

            return None

        print(
            "✅ Player count is correct."
        )

        # ====================================================
        # SEPARATE HIGH-LEVEL AND REGULAR PLAYERS
        # ====================================================

        high_level_players, regular_players = (
            separate_players(all_players)
        )

        print()
        print("=" * 80)
        print("AUTOMATIC PLAYER CLASSIFICATION")
        print("=" * 80)

        print(
            f"Total players      : {len(all_players)}"
        )

        print(
            f"High-level players : "
            f"{len(high_level_players)}"
        )

        print(
            f"Regular players    : "
            f"{len(regular_players)}"
        )

        # ====================================================
        # HIGH-LEVEL PLAYERS BY POSITION
        # ====================================================

        high_level_by_position = group_by_position(
            high_level_players
        )

        print()
        print("HIGH-LEVEL PLAYERS BY POSITION")

        for position in POSITIONS:

            print(
                f"{position:<15}: "
                f"{len(high_level_by_position[position])}"
            )

        # ====================================================
        # STEP 1:
        # DISTRIBUTE HIGH-LEVEL PLAYERS BY POSITION
        # ====================================================

        distribute_high_level_players(
            high_level_players,
            Black_team,
            White_team
        )

        # ====================================================
        # STEP 2:
        # DISTRIBUTE REGULAR PLAYERS
        # ====================================================

        distribute_regular_players(
            regular_players,
            Black_team,
            White_team,
            all_players
        )

        # ====================================================
        # FINAL RANDOMIZATION
        # ====================================================

        random.shuffle(Black_team)
        random.shuffle(White_team)

        # ====================================================
        # DISPLAY RESULTS
        # ====================================================

        display_high_level_summary(
            all_players,
            Black_team,
            White_team
        )

        display_high_level_detail(
            Black_team,
            White_team
        )

        display_high_level_position_balance(
            Black_team,
            White_team
        )

        display_team(
            Black_team,
            "BLACK TEAM"
        )

        display_team(
            White_team,
            "WHITE TEAM"
        )

        # ====================================================
        # TEAM TOTALS
        # ====================================================

        print()
        print("=" * 80)
        print("TEAM TOTALS")
        print("=" * 80)

        print(
            f"Black Team : {len(Black_team)} players"
        )

        print(
            f"White Team : {len(White_team)} players"
        )

        # ====================================================
        # POSITION SUMMARY
        # ====================================================

        display_position_summary(
            Black_team,
            White_team
        )

        # ====================================================
        # BALANCE CHECK
        # ====================================================

        display_balance_check(
            Black_team,
            White_team
        )

        # ====================================================
        # RETURN TEAMS
        # ====================================================

        return Black_team, White_team

    except KeyboardInterrupt:

        print()

        print(
            "🛑 Program cancelled by user."
        )

        return None

    except Exception as error:

        print()
        print("=" * 80)
        print("UNEXPECTED ERROR")
        print("=" * 80)

        print(
            f"Error type: {type(error).__name__}"
        )

        print(
            f"Error: {error}"
        )

        print(
            "The program stopped safely."
        )

        return None


# ============================================================
# PROGRAM ENTRY POINT
# ============================================================

if __name__ == "__main__":

    Players_Grouping()